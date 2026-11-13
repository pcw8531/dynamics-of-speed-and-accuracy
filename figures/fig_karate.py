"""Prediction 2 in Kyokushin karate, the test as run (analysis/plan_karate.md, analysis/karate_test.py).

A, each athlete's H (degree shares of the 42 strongest links among 14 joints, mean of S01 and S02) against grade.
B, the primary test: rank of H and rank of grade, each with the rank of age removed, and their correlation.
Reads analysis/karate_per_athlete.csv. Writes figures/fig_karate.png.
"""
import logging
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
ROOT = Path(__file__).resolve().parent.parent
mpl.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 7.5,
                     "mathtext.fontset": "stix", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 1.8, "ytick.major.size": 1.8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 7.5, "savefig.dpi": 600,
                     "axes.spines.top": False, "axes.spines.right": False})
INK, GREY = "#000000", "#6e6e6e"

ath = pd.read_csv(ROOT / "analysis" / "karate_per_athlete.csv")
g, H, age = ath.grade_rank.values, ath.H.values, ath.age_years.values

fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.3))
fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.2, wspace=0.32)

sc = a.scatter(g + np.linspace(-0.12, 0.12, len(g))[np.argsort(np.argsort(age))], H, s=16, c=age, cmap="Greys",
               vmin=0, vmax=55, edgecolors=INK, linewidths=0.5, zorder=3)
a.set_xticks([1, 5, 9, 13]); a.set_xticklabels(["9 kyu", "5 kyu", "1 kyu", "4 dan"])
a.set_xlim(0.3, 13.7); a.set_xlabel("grade"); a.set_ylabel(r"$H$, bits")
cb = fig.colorbar(sc, ax=a, fraction=0.05, pad=0.02); cb.set_label("age, years", fontsize=7); cb.ax.tick_params(labelsize=6.5)
cb.outline.set_linewidth(0.5)
a.set_title("A", loc="left", fontsize=9, fontweight="bold")

ra, rg, rh = rankdata(age), rankdata(g), rankdata(H)
X = np.c_[np.ones_like(ra), ra]
eg = rg - X @ np.linalg.lstsq(X, rg, rcond=None)[0]; eh = rh - X @ np.linalg.lstsq(X, rh, rcond=None)[0]
r = np.corrcoef(eg, eh)[0, 1]
b.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2))); b.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
b.scatter(eg, eh, s=16, color="white", edgecolors=INK, linewidths=0.6, zorder=3)
xx = np.array([eg.min(), eg.max()]); b.plot(xx, np.polyfit(eg, eh, 1)[0] * xx, color=INK, lw=1.0)
b.set_title(f"$r$ = {r:.3f}, two-sided $p$ = .048", loc="right", fontsize=7)
b.set_xlabel("rank of grade, age removed"); b.set_ylabel(r"rank of $H$, age removed")
b.set_title("B", loc="left", fontsize=9, fontweight="bold")
fig.savefig(ROOT / "figures" / "fig_karate.png")
print(f"r = {r:.4f} from {len(ath)} athletes")
