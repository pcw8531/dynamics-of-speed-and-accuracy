"""Prediction 2 test in Kyokushin karate, run once, as fixed in analysis/plan_karate.md (with its recorded changes).

Usage: python -I analysis/karate_test.py <folder with the athlete folders of Karate.zip> <clone of py-c3d at b3060462bf35>
Reads S01 and S02, condition E01, trials T01 and T02 of every athlete. Writes karate_per_recording.csv and
karate_per_athlete.csv next to this script and prints the test statistics.
"""
import re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
from scipy.stats import rankdata, spearmanr

DATA, PYC3D = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(PYC3D))
import c3d                                                            # noqa: E402  py-c3d, MIT licence
warnings.filterwarnings("ignore", message="No analog data found in file.")

HERE = Path(__file__).resolve().parent
TABLE1 = HERE.parent / "data" / "transcribed" / "szczesna2021_table1_participants.csv"
JOINTS = ["LHipAngles", "RHipAngles", "LKneeAngles", "RKneeAngles", "LAnkleAngles", "RAnkleAngles",
          "LShoulderAngles", "RShoulderAngles", "LElbowAngles", "RElbowAngles", "LWristAngles", "RWristAngles",
          "LSpineAngles", "LNeckAngles"]                             # rule 1
N = len(JOINTS); LINKS = 3 * N; CUT_HZ = 15.0; MIN_VALID = 0.80; MAX_GAP_S = 0.1
SEED, PERMS = 20261105, 10_000
PAT = re.compile(r"^(\d{4}-\d{2}-\d{2})-(B\d{4})-(S01|S02)-E01-(T01|T02)\.c3d$")


def read(path):
    """Frame rate from the header and angles (frames, 14, 3) in degrees. A frame is invalid for a joint when the point
    is flagged missing (residual below 0) or its three components are all exactly zero (rule 2)."""
    with open(path, "rb") as fh:
        r = c3d.Reader(fh)
        labels = [s.strip() for s in r.point_labels]
        idx = [labels.index(j) for j in JOINTS]                       # every output is present (recorded change)
        P = np.array([pts[:, :4].copy() for _, pts, _ in r.read_frames()])
        rate = float(r.header.frame_rate)
    A = P[:, idx, :3].copy()
    bad = (P[:, idx, 3] < 0) | np.all(A == 0, axis=2)
    A[bad] = np.nan
    return rate, A


def clean(A, rate):
    """Rule 2 with the recorded change: at least 80 percent of frames valid for all joints, edge frames trimmed,
    gaps inside filled linearly up to 0.1 s, a recording with a longer gap not used."""
    ok = ~np.isnan(A).any(axis=(1, 2))
    info = dict(frames=len(ok), valid=float(ok.mean()))
    if info["valid"] < MIN_VALID:
        return None, info | dict(reason="under 80 percent valid")
    first, last = np.flatnonzero(ok)[[0, -1]]
    A, ok = A[first:last + 1], ok[first:last + 1]
    info["trimmed"] = int(info["frames"] - len(ok))
    gaps, i = [], 0
    while i < len(ok):
        if not ok[i]:
            j = i
            while not ok[j]:
                j += 1
            gaps.append(j - i); i = j
        else:
            i += 1
    info["max_gap"] = max(gaps, default=0)
    if info["max_gap"] > round(MAX_GAP_S * rate):
        return None, info | dict(reason="gap over 0.1 s")
    if gaps:
        t = np.arange(len(ok))
        for k in range(N):
            for c in range(3):
                A[~ok, k, c] = np.interp(t[~ok], t[ok], A[ok, k, c])
    return A, info


def speeds(A, rate):
    """Rule 3: unwrap in degrees, fourth-order zero-lag Butterworth low-pass at 15 Hz, angular speed in degrees per
    second as the norm of the first difference of the three angles times the rate."""
    b, a = butter(4, CUT_HZ / (rate / 2))
    U = np.unwrap(A, period=360.0, axis=0)
    F = filtfilt(b, a, U, axis=0)
    return np.linalg.norm(np.diff(F, axis=0), axis=2).T * rate      # (14, frames - 1)


def entropy(p):
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def coupling_H(S):
    """Rule 4: absolute Pearson correlations, the 3N strongest of the 91 pairs (ties by a stable sort of pair index),
    H of the degree shares. Also the strength-share H of the full network (secondary c)."""
    W = np.abs(np.corrcoef(S))
    iu = np.triu_indices(N, 1)
    w = W[iu]
    keep = np.argsort(-w, kind="stable")[:LINKS]
    k = np.zeros(N)
    np.add.at(k, iu[0][keep], 1); np.add.at(k, iu[1][keep], 1)
    np.fill_diagonal(W, 0.0)
    s = W.sum(1)
    return entropy(k / k.sum()), entropy(s / s.sum()), int((k == 0).sum())


def partial_spearman(grade, H, age, seed=SEED):
    """Primary test: ranks of grade and H each regressed on the rank of age, r of the residuals, one-sided p from
    10,000 permutations of the H residuals (and two-sided p for the rule against)."""
    rg, rh, ra = rankdata(grade), rankdata(H), rankdata(age)
    X = np.c_[np.ones_like(ra), ra]
    eg = rg - X @ np.linalg.lstsq(X, rg, rcond=None)[0]
    eh = rh - X @ np.linalg.lstsq(X, rh, rcond=None)[0]
    r = float(np.corrcoef(eg, eh)[0, 1])
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(eg, rng.permutation(eh))[0, 1] for _ in range(PERMS)])
    p_one = (1 + np.sum(null <= r)) / (PERMS + 1)
    p_two = (1 + np.sum(np.abs(null) >= abs(r))) / (PERMS + 1)
    verdict = "support" if (r < 0 and p_one < .05) else "against" if (r > 0 and p_two < .05) else "undecided"
    return r, float(p_one), float(p_two), verdict


def main():
    t1 = pd.read_csv(TABLE1).set_index("code")
    rows = []
    for f in sorted(DATA.glob("B*/*/*.c3d")):
        m = PAT.match(f.name)
        if not m:
            continue
        rate, A = read(f)
        A, info = clean(A, rate)
        row = dict(file=f.name, athlete=m.group(2), technique=m.group(3), trial=m.group(4), rate=rate, **info)
        if A is not None:
            row["H_degree"], row["H_strength"], row["isolated_joints"] = coupling_H(speeds(A, rate))
            row["used"] = True
        else:
            row["used"] = False
        rows.append(row)
    rec = pd.DataFrame(rows)
    rec.to_csv(HERE / "karate_per_recording.csv", index=False, float_format="%.6f")

    use = rec[rec.used]
    tech = use.groupby(["athlete", "technique"])[["H_degree", "H_strength"]].mean().unstack("technique")
    ath = pd.DataFrame({"H_S01": tech[("H_degree", "S01")], "H_S02": tech[("H_degree", "S02")],
                        "Hs_S01": tech[("H_strength", "S01")], "Hs_S02": tech[("H_strength", "S02")]})
    ath["H"] = ath[["H_S01", "H_S02"]].mean(axis=1)                   # rule 5
    ath["H_strength"] = ath[["Hs_S01", "Hs_S02"]].mean(axis=1)
    ath["recordings"] = use.groupby("athlete").size()
    ath = ath.join(t1[["age_years", "grade", "grade_rank", "training_years"]])
    ath.index.name = "athlete"
    ath.to_csv(HERE / "karate_per_athlete.csv", float_format="%.6f")
    left_out = sorted(set(t1.index) - set(ath.index))

    print("Prediction 2, Kyokushin karate, joint coupling (analysis/plan_karate.md)")
    print(f"recordings read {len(rec)}, used {int(rec.used.sum())}, not used {int((~rec.used).sum())}"
          + ("" if rec.used.all() else f" {rec.loc[~rec.used, ['file', 'reason']].values.tolist()}"))
    print(f"rates {rec.rate.value_counts().to_dict()}; edge frames trimmed in {int((rec.get('trimmed', 0) > 0).sum())} recordings; "
          f"interpolated gaps in {int((rec.get('max_gap', 0) > 0).sum())}")
    print(f"athletes {len(ath)} of {len(t1)}, left out {left_out}; athletes with one technique only "
          f"{int(ath[['H_S01', 'H_S02']].isna().any(axis=1).sum())}")
    print(f"recordings with a joint outside the 42 strongest links: {int((use.isolated_joints > 0).sum())}")
    print(f"H, degree shares of the 42 strongest links: mean {ath.H.mean():.4f}, sd {ath.H.std(ddof=1):.4f}, "
          f"range {ath.H.min():.4f} to {ath.H.max():.4f}; ceiling log2 14 = {np.log2(N):.4f}")
    r, p1, p2, v = partial_spearman(ath.grade_rank.values, ath.H.values, ath.age_years.values)
    print(f"PRIMARY  partial Spearman of grade rank and H, age controlled: r = {r:.4f}, one-sided p = {p1:.4f}, "
          f"two-sided p = {p2:.4f}  ->  {v.upper()}")
    rho, pa = spearmanr(ath.grade_rank, ath.H)
    print(f"secondary (a) Spearman of grade rank and H without age: rho = {rho:.4f}, two-sided p = {pa:.4f}")
    for t in ("S01", "S02"):
        sub = ath.dropna(subset=[f"H_{t}"])
        r_, p1_, p2_, v_ = partial_spearman(sub.grade_rank.values, sub[f"H_{t}"].values, sub.age_years.values)
        print(f"secondary (b) {t} only, n = {len(sub)}: r = {r_:.4f}, one-sided p = {p1_:.4f}, two-sided p = {p2_:.4f} ({v_})")
    r_, p1_, p2_, v_ = partial_spearman(ath.grade_rank.values, ath.H_strength.values, ath.age_years.values)
    print(f"secondary (c) strength shares of the full network: mean H {ath.H_strength.mean():.4f}, "
          f"r = {r_:.4f}, one-sided p = {p1_:.4f}, two-sided p = {p2_:.4f} ({v_})")
    print(f"grade and age: Spearman {spearmanr(ath.grade_rank, ath.age_years)[0]:.3f}")


if __name__ == "__main__":
    main()
