"""Prediction 4 test in How We Type, as fixed in analysis/plan_P4_howwetype.md.

Usage: python -I analysis/p4_howwetype.py <typing folder of *_matched.txt>
Reads H and the median inter-key interval from analysis/p2_per_typist.csv, written by p2_test.py.
"""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

SEED, N_BOOT = 20261104, 10_000
here = Path(__file__).resolve().parent
per = pd.read_csv(here / "p2_per_typist.csv")

corr = {}
for f in sorted(glob.glob(str(Path(sys.argv[1]) / "*_matched.txt"))):
    t = pd.read_csv(f, sep="\t")
    n_bs = int((t.key_symbol == "BackSpace").sum())
    final = t.groupby("stimulus", sort=False).transcribed.last().astype(str)
    corr[int(t.user_id.iloc[0])] = n_bs / int(final.str.len().sum())
per["corrections_per_char"] = per.uid.map(corr)

rng = np.random.default_rng(SEED)
out = {}
for name in ("median_iki_ms", "corrections_per_char"):
    r, p2 = spearmanr(per.H_bits, per[name])
    p1 = p2 / 2 if r > 0 else 1 - p2 / 2
    boot = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(per), len(per))
        boot.append(spearmanr(per.H_bits.values[i], per[name].values[i])[0])
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    out[name] = (r, p1, p2, lo, hi)
    print(f"rho(H, {name}) = {r:.3f}, one-sided p = {p1:.4f}, two-sided p = {p2:.4f}, 95% interval [{lo:.3f}, {hi:.3f}]")
sup = all(r > 0 and p1 < .05 for r, p1, p2, lo, hi in out.values())
ag = any(r < 0 and p2 < .05 for r, p1, p2, lo, hi in out.values())
print("verdict:", "support" if sup else ("against" if ag else "undecided"))
per.to_csv(here / "p4_per_typist.csv", index=False)
