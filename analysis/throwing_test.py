"""Throwing practice test, run once, as fixed in analysis/plan_throwing.md (with its recorded change).

Usage: python -I analysis/throwing_test.py <folder with S1.zip ... S7.zip and Trialinfo.xlsx>
Reads the trial CSVs directly from the zips. Writes throwing_per_throw.csv and throwing_per_session.csv
next to this script and prints the test statistics.
"""
import io, re, sys, zipfile, collections
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
from scipy.stats import wilcoxon, spearmanr

D = Path(sys.argv[1]); FS = 200.0
NONBODY = {"Ball", "FP1", "FP2"}
PAT = re.compile(r"S(\d)_(DP|NDP)_(\d+)_+([A-Z]+)\.csv$")
SESSIONS = [f"{w}w{d}d" for w in range(1, 6) for d in range(1, 4)]
b, a = butter(4, 15.0 / (FS / 2))


def ball_speeds():
    out = {}
    x = pd.ExcelFile(D / "Trialinfo.xlsx")
    for s in range(1, 8):
        df = x.parse(f"S{s}", header=None)
        for c in range(0, df.shape[1], 3):
            sess = None
            for r in range(df.shape[0]):
                v = df.iat[r, c]
                if isinstance(v, str) and re.fullmatch(r"\dw\dd", v.strip()):
                    sess = v.strip(); continue
                if isinstance(v, str) and sess and re.fullmatch(r"N?DP_\d+", v.strip()):
                    side, k = v.strip().split("_")
                    sp = pd.to_numeric(df.iat[r, c + 2], errors="coerce")
                    out[(s, sess, side, int(k))] = float(sp) if pd.notna(sp) else np.nan
    return out


def read_trial(z, name):
    lines = z.read(name).decode("latin-1").splitlines()
    lab = lines[1].split("\t")
    markers = [lab[i] for i in range(4, len(lab), 3) if lab[i]]
    df = pd.read_csv(io.StringIO("\n".join(lines[5:])), sep="\t", header=None)
    ev = df.iloc[0, 1:4].astype(float).values
    if np.isnan(ev).any() or not (ev[0] < ev[1] < ev[2]):
        return None
    fr = df.iloc[:, 0].astype(float).values
    w = (fr >= round(ev[0] * FS) + 1) & (fr <= round(ev[2] * FS) + 1)
    X = df.iloc[:, 4:4 + 3 * len(markers)].apply(pd.to_numeric, errors="coerce").values[w]
    return {mk: X[:, 3 * j:3 * j + 3] for j, mk in enumerate(markers) if mk not in NONBODY}


def entropy(P):
    """P: dict marker -> positions (frames x 3) inside the window, all complete."""
    names = sorted(P)
    S = []
    for mk in names:
        Q = filtfilt(b, a, P[mk], axis=0)
        S.append(np.linalg.norm(np.diff(Q, axis=0), axis=1) * FS)
    S = np.array(S)
    W = np.abs(np.corrcoef(S)); np.fill_diagonal(W, 0.0)
    N = len(names); iu = np.triu_indices(N, 1)
    order = np.argsort(W[iu])[::-1][:3 * N]
    k = np.zeros(N)
    for t in order:
        k[iu[0][t]] += 1; k[iu[1][t]] += 1
    p = k[k > 0] / k.sum()
    return float(-(p * np.log2(p)).sum()), N


# pass 1: which markers are present, and their gap counts, over the clean throws
trials = []
for s in range(1, 8):
    z = zipfile.ZipFile(D / f"S{s}.zip")
    for name in z.namelist():
        m = PAT.search(name)
        if m and m.group(4) == "O":
            trials.append((s, name.split("/")[-2], m.group(2), int(m.group(3)), name))
present, gaps, cache = [], collections.Counter(), {}
for s, sess, side, k, name in trials:
    P = read_trial(zipfile.ZipFile(D / f"S{s}.zip"), name)
    if P is None:
        continue
    present.append(set(P))
    for mk, x in P.items():
        if np.isnan(x).any():
            gaps[mk] += 1
    cache[(s, sess, side, k)] = name
n = len(present)
body = set.intersection(*present)
keep = sorted(mk for mk in body if gaps[mk] <= 0.02 * n)
print(f"clean throws readable {n}; body markers in all {len(body)}; kept {len(keep)}; dropped {sorted(body - set(keep))}")

# pass 2: H for every throw that has no gap in the kept markers
speeds = ball_speeds(); rows = []
for (s, sess, side, k), name in cache.items():
    P = read_trial(zipfile.ZipFile(D / f"S{s}.zip"), name)
    P = {mk: P[mk] for mk in keep}
    if any(np.isnan(x).any() for x in P.values()):
        continue
    H, N = entropy(P)
    rows.append((s, sess, SESSIONS.index(sess) + 1, side, k, H, speeds.get((s, sess, side, k), np.nan)))
t = pd.DataFrame(rows, columns=["subject", "session", "session_no", "side", "trial", "H_bits", "ball_speed"])
here = Path(__file__).resolve().parent
t.to_csv(here / "throwing_per_throw.csv", index=False)
print(f"throws analysed {len(t)}; N = {len(keep)} markers; links per network {3 * len(keep)}; ceiling log2 N = {np.log2(len(keep)):.4f}")

ses = t.groupby(["subject", "side", "session_no"]).agg(H=("H_bits", "mean"), speed=("ball_speed", "mean"), n=("H_bits", "size")).reset_index()
ses.to_csv(here / "throwing_per_session.csv", index=False)

slope = lambda g, col: np.polyfit(g.session_no, g[col], 1)[0]
sl = ses.groupby(["subject", "side"]).apply(lambda g: slope(g, "H")).unstack()
sl_sp = ses.groupby(["subject", "side"]).apply(lambda g: slope(g.dropna(subset=["speed"]), "speed")).unstack()
per_subject = sl.mean(axis=1)
w1 = wilcoxon(per_subject, alternative="less"); w1b = wilcoxon(per_subject)
print("\nPrimary: slope of session mean H per session (bits), per subject and side")
print(sl.round(5).to_string())
print("mean over sides:", per_subject.round(5).tolist())
print(f"one-sided Wilcoxon (slope < 0): W = {w1.statistic}, p = {w1.pvalue:.4f}; two-sided p = {w1b.pvalue:.4f}")
verdict = "support" if w1.pvalue < .05 else ("against" if (w1b.pvalue < .05 and per_subject.median() > 0) else "undecided")
print("verdict:", verdict)

rs = t.dropna(subset=["ball_speed"]).groupby(["subject", "side"]).apply(lambda g: spearmanr(g.H_bits, g.ball_speed)[0]).unstack()
rsub = rs.mean(axis=1); w2 = wilcoxon(rsub, alternative="less")
print("\nSecondary: Spearman rho(H, ball speed) per subject and side"); print(rs.round(3).to_string())
print(f"mean over sides: {rsub.round(3).tolist()}; one-sided Wilcoxon (rho < 0): p = {w2.pvalue:.4f}")

mh = ses.groupby(["subject", "side"]).H.mean().unstack()
w3 = wilcoxon(mh["DP"], mh["NDP"], alternative="less")
print("\nSecondary: mean H over sessions, dominant against non-dominant"); print(mh.round(4).to_string())
print(f"one-sided paired Wilcoxon (DP < NDP): p = {w3.pvalue:.4f}")

print("\nPractice check: slope of session mean ball speed per session"); print(sl_sp.round(3).to_string())
