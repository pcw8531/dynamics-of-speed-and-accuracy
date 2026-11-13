"""Prediction 5 test in the Axon game data, as fixed in analysis/plan_P5_axon.md.
Usage: python -I analysis/p5_axon.py <data_by_cookie.json>"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest
from scipy.optimize import curve_fit

d = json.load(open(sys.argv[1]))
S = np.array([[v[k][0] for k in sorted(v)[:30]] for v in d.values() if 30 <= len(v) <= 100], float)
m = S.reshape(len(S), 3, 10).mean(axis=2)
g1, g2 = m[:, 1] - m[:, 0], m[:, 2] - m[:, 1]
dec = g1 - g2
k, n = int((dec > 0).sum()), len(dec)
bt = binomtest(k, n, 0.5, alternative="greater"); bt2 = binomtest(k, n, 0.5)
lo, hi = bt2.proportion_ci(confidence_level=0.95, method="wilson")
share = k / n
if share > 0.5 and bt.pvalue < .05 and np.median(g1) > 0:
    verdict = "support"
elif (share < 0.5 and bt2.pvalue < .05) or np.median(g1) <= 0:
    verdict = "against"
else:
    verdict = "undecided"
print(f"identities {n}; share with d > 0 = {share:.4f} [{lo:.4f}, {hi:.4f}], one-sided p = {bt.pvalue:.3g}")
print(f"median g1 = {np.median(g1):.1f}, median g2 = {np.median(g2):.1f}, median d = {np.median(dec):.1f}")
print("verdict:", verdict)

f = lambda t, Sinf, R, a: Sinf - R * (1 + t) ** (-a)
t = np.arange(30)
def fit(M):
    return curve_fit(f, t, M, p0=[M[-1], M[-1] - M[0], 0.5], maxfev=20000)[0]
p = fit(S.mean(axis=0))
rng = np.random.default_rng(20261104)
boot = []
for _ in range(1000):
    try:
        boot.append(fit(S[rng.integers(0, n, n)].mean(axis=0))[2])
    except RuntimeError:
        pass
a_lo, a_hi = np.percentile(boot, [2.5, 97.5])
print(f"average curve, Equation 18: S_inf = {p[0]:.0f}, R = {p[1]:.0f}, alpha = {p[2]:.3f} [{a_lo:.3f}, {a_hi:.3f}] ({len(boot)} resamples converged)")
np.savetxt(Path(__file__).with_name("p5_mean_curve.csv"), np.c_[t + 1, S.mean(axis=0), f(t, *p)], delimiter=",",
           header="play,mean_score,fit_eq18", comments="", fmt=["%d", "%.2f", "%.2f"])
