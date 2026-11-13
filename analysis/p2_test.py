"""Prediction 2 test, run once on the full data, as fixed in analysis/plan_P2.md.

Usage:
  python -I analysis/p2_test.py <typing folder of *_matched.txt> <keyboard_flat_coordinates.csv> <motion capture folder>

The reader for the motion capture files (read_mocap) follows the column names in the dataset readme. If the
released files differ in layout, only the reader is adapted; the rules in FIXED stay as they are.
"""
import glob, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import butter, filtfilt

FIXED = dict(letters_only=True, max_iki_ms=1000, key_width_mm=19.0, gap_fill_frames=24, min_stretch_frames=240,
             lowpass_hz=10.0, fs_hz=240.0, marker_valid_share=0.80, n_boot=10_000, seed=20261104,
             equivalence=(0.80, 1.25))
LETTERS = set("abcdefghijklmnopqrstuvwxyz") | {"aring", "adiaeresis", "odiaeresis"}


def moves(typing_dir, coord_csv):
    c = pd.read_csv(coord_csv, sep="\t")
    xyz = 1000 * c.drop_duplicates("key").set_index("key")[["x", "y", "z"]]
    rows, speed = [], {}
    for f in sorted(glob.glob(str(Path(typing_dir) / "*_matched.txt"))):
        t = pd.read_csv(f, sep="\t")
        uid = int(t.user_id.iloc[0])
        k, fi, iki, st = t.key_symbol.values, t.finger.values, t.iki.values.astype(float), t.stimulus.values
        first = np.r_[True, st[1:] != st[:-1]]
        keep = ~first & (iki > 0) & (iki <= FIXED["max_iki_ms"])
        speed[uid] = float(np.median(iki[keep]))
        bs = k == "BackSpace"
        near_bs = np.r_[False, bs[:-1]] | np.r_[bs[1:], False]          # BackSpace just before or just after
        for i in np.where(keep)[0]:
            a, b = k[i - 1], k[i]
            if a == b or fi[i] != fi[i - 1] or a not in LETTERS or b not in LETTERS:
                continue
            if near_bs[i] or near_bs[i - 1]:
                continue
            D = float(np.linalg.norm(xyz.loc[b].values - xyz.loc[a].values))
            rows.append((uid, np.log2(2 * D / FIXED["key_width_mm"]), iki[i]))
    return pd.DataFrame(rows, columns=["uid", "ID", "iki"]), pd.Series(speed)


def read_mocap(path):
    """The released files are tab-separated, with one line of capture settings and one blank line before
    the header. Quote marks inside the stimulus text are literal characters, not field delimiters."""
    m = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False, quoting=3)
    m = m[m.condition == "sentences"]
    cols = [c for c in m.columns if re.match(r"^Hands_[LR]_.+_[xyz]$", c)]
    return int(m.uid.iloc[0]), m, cols


def marker_speeds(m, cols):
    """Speed series of each marker within typing frames, by the gap, filter and stretch rules in FIXED."""
    b, a = butter(4, FIXED["lowpass_hz"] / (FIXED["fs_hz"] / 2))
    typing = np.zeros(len(m), bool)                                  # m has a 0..n-1 index
    for _, g in m.groupby("stimulus_index"):
        k = g.index[g.input_index.notna().values]
        if len(k):
            typing[k[0]:k[-1] + 1] = True
    names = sorted({c[:-2] for c in cols})
    out = {}
    for n in names:
        P = m[[n + "_x", n + "_y", n + "_z"]].interpolate(limit=FIXED["gap_fill_frames"], limit_area="inside").values
        ok = ~np.isnan(P).any(axis=1)
        V = np.full(len(P), np.nan)
        edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(int), 0]))
        for s, e in zip(edges[::2], edges[1::2]):
            if e - s >= FIXED["min_stretch_frames"]:
                Q = filtfilt(b, a, P[s:e], axis=0)
                V[s + 1:e] = np.linalg.norm(np.diff(Q, axis=0), axis=1) * FIXED["fs_hz"]
        out[n] = V[typing]
    return pd.DataFrame(out)


def entropy_of(S):
    W = S.corr(min_periods=int(FIXED["fs_hz"])).abs().to_numpy(copy=True)
    np.fill_diagonal(W, 0.0)
    s = W.sum(axis=1)
    p = s / s.sum()
    return float(-(p * np.log2(p)).sum())


def group_slope(by_uid, uids):
    """One slope for the listed typists, one intercept for each entry. A typist drawn twice in a
    bootstrap resample enters twice, with its own intercept each time."""
    ID = np.concatenate([by_uid[u].ID.values for u in uids])
    y = np.concatenate([by_uid[u].iki.values for u in uids])
    lab = np.concatenate([np.full(len(by_uid[u]), j) for j, u in enumerate(uids)])
    X = np.column_stack([ID] + [(lab == j).astype(float) for j in range(len(uids))])
    return float(np.linalg.lstsq(X, y, rcond=None)[0][0])


def main(typing_dir, coord_csv, mocap_dir):
    mv, speed = moves(typing_dir, coord_csv)
    by_uid = {u: g for u, g in mv.groupby("uid")}
    speeds = {}
    for f in sorted(glob.glob(str(Path(mocap_dir) / "*_mocap_matched.csv"))):
        uid, m, cols = read_mocap(f)
        speeds[uid] = marker_speeds(m.reset_index(drop=True), cols)
    valid = {u: S.notna().mean() for u, S in speeds.items()}
    markers = sorted(set.intersection(*[set(v.index[v >= FIXED["marker_valid_share"]]) for v in valid.values()]))
    H = pd.Series({u: entropy_of(S[markers]) for u, S in speeds.items()})
    uids = sorted(set(speed.index) & set(H.index))
    order = speed[uids].sort_values().index
    fast, slow = list(order[:len(order) // 2]), list(order[len(order) // 2:])

    def stat(f_ids, s_ids):
        bf, bs = group_slope(by_uid, f_ids), group_slope(by_uid, s_ids)
        hf, hs = H[f_ids].mean(), H[s_ids].mean()
        return bf, bs, hf, hs, (bf / bs) / (hf / hs)

    bf, bs, hf, hs, rho = stat(fast, slow)
    rng = np.random.default_rng(FIXED["seed"])
    boot = np.array([stat(list(rng.choice(fast, len(fast))), list(rng.choice(slow, len(slow))))
                     for _ in range(FIXED["n_boot"])])
    ci = lambda x: np.percentile(x, [2.5, 97.5])
    lo, hi = ci(boot[:, 4]); bf_ci, bs_ci = ci(boot[:, 0]), ci(boot[:, 1])
    lower, upper = FIXED["equivalence"]
    if bf_ci[0] <= 0 <= bf_ci[1] or bs_ci[0] <= 0 <= bs_ci[1]:
        verdict = "undecided (a group slope interval includes 0)"
    elif lo > 1 or hi < 1:
        verdict = "against (interval excludes 1)"
    elif lower <= lo and hi <= upper:
        verdict = "support (interval includes 1 and lies within 0.80 to 1.25)"
    else:
        verdict = "undecided (interval includes 1 but is wider than 0.80 to 1.25)"
    print(f"moves kept {len(mv)}, typists {len(uids)}, markers {len(markers)}")
    print(f"faster: b = {bf:.2f} ms/bit [{bf_ci[0]:.2f}, {bf_ci[1]:.2f}], H = {hf:.4f} bits, beta = {bf / hf:.3f}")
    print(f"slower: b = {bs:.2f} ms/bit [{bs_ci[0]:.2f}, {bs_ci[1]:.2f}], H = {hs:.4f} bits, beta = {bs / hs:.3f}")
    print(f"rho = {rho:.3f}, 95% interval [{lo:.3f}, {hi:.3f}]: {verdict}")
    per = pd.DataFrame({"uid": uids, "median_iki_ms": speed[uids].values, "H_bits": H[uids].values,
                        "slope_ms_per_bit": [group_slope(by_uid, [u]) for u in uids],
                        "group": ["faster" if u in fast else "slower" for u in uids]})
    per.to_csv(Path(__file__).with_name("p2_per_typist.csv"), index=False)


if __name__ == "__main__":
    main(*sys.argv[1:4])
