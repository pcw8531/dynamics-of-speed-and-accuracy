"""Figure S4 of the Supplemental Material: skill and practice in two open records, from the analyses run under the plans of the
analysis folder (Supplemental Section S7).

A, How We Type (Feit et al., 2016a, 2016b): the slope of the inter-key interval on the index of difficulty for
same-finger moves, with one intercept for each typist and one slope for each half of the 30 typists split by median
inter-key interval, with the 95 percent bootstrap intervals of the plan (analysis/p2_result.txt, analysis/p2_per_typist.csv).
B, the online game of Stafford and Dewar (2014a, 2014b): mean score over the first 30 plays of the 5,293 players with 30
to 100 plays, with the least-squares fit of Equation 16 written for a score that rises, S(t) = S_inf - R (1 + t)^-alpha
(analysis/p5_mean_curve.csv, analysis/p5_axon_result.txt). Palette and type from figures/style.py.
Run from the repository root: python figures/figS4_open_records.py
"""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, EDGE, RULE, GREY, TONE, tb_axes, save_panels  # noqa: E402

A = ROOT / "analysis"
txt = (A / "p2_result.txt").read_text()
grp = {}
for g in ("faster", "slower"):
    m = re.search(g + r": b = ([\d.]+) ms/bit \[([\d.]+), ([\d.]+)\]", txt)
    grp[g] = tuple(float(v) for v in m.groups())
typ = pd.read_csv(A / "p2_per_typist.csv")
iki = typ.groupby("group").median_iki_ms.median()
curve = pd.read_csv(A / "p5_mean_curve.csv")
ax5 = (A / "p5_axon_result.txt").read_text()
n_id = int(re.search(r"identities (\d+)", ax5).group(1))
share = float(re.search(r"share with d > 0 = ([\d.]+)", ax5).group(1))
alpha = re.search(r"alpha = ([\d.]+) \[([\d.]+), ([\d.]+)\]", ax5).groups()

S = Sheet(6.5, 1.95)
axA = S.ax(0.62, 0.40, 2.25, 1.30)
axB = S.ax(3.98, 0.40, 2.25, 1.30)

# ---------------------------------------------------------------- A, the slope of the faster and the slower typists
ax = axA
for x_, g in ((0, "slower"), (1, "faster")):
    b, lo, hi = grp[g]
    ax.plot([x_, x_], [lo, hi], color=EDGE, lw=0.9, zorder=2)
    ax.plot(x_, b, "o", ms=5.0, mfc=TONE["random"] if g == "slower" else INK,
            mec="white", mew=0.5, zorder=3)
    ax.text(x_ + 0.08, b, f"{b:.0f}", fontsize=6.6, ha="left", va="center")
ax.annotate("", xy=(0.92, grp["faster"][0]), xytext=(0.08, grp["slower"][0]),
            arrowprops=dict(arrowstyle="-|>", mutation_scale=7, lw=0.7, color=GREY, ls=(0, (2.0, 1.5))))
ax.text(0.5, (grp["faster"][0] + grp["slower"][0]) / 2 + 10, f"ratio {grp['faster'][0] / grp['slower'][0]:.2f}",
        fontsize=6.4, ha="center", va="bottom", color=EDGE)
ax.set_xlim(-0.45, 1.45); ax.set_xticks([0, 1])
ax.set_xticklabels([f"slower half\n{iki['slower']:.0f} ms", f"faster half\n{iki['faster']:.0f} ms"], fontsize=6.6)
tb_axes(ax, yticks=[0, 40, 80], ylim=(0, 84), grid="y")
ax.set_ylabel("slope, ms per bit", labelpad=1.5)
ax.text(-0.40, 6, "15 typists each, median\ninter-key interval below", fontsize=5.9, ha="left", va="bottom", color=GREY)

# ---------------------------------------------------------------- B, the course of practice in the game
ax = axB
t = curve.play.values
ax.plot(t, curve.fit_eq18 / 1000, color=INK, lw=1.0, zorder=2)
ax.plot(t, curve.mean_score / 1000, "o", ls="", ms=3.0, mfc=TONE["small-world"], mec=EDGE, mew=0.45, zorder=3)
ax.text(29.5, 14.5, f"{n_id:,} players\n{100 * share:.1f} percent gained more\nearly than later\n"
        + rf"$\alpha$ = {alpha[0]} ({alpha[1]} to {alpha[2]})", fontsize=6.2, ha="right", va="bottom", linespacing=1.15)
tb_axes(ax, xticks=[1, 15, 30], xlim=(0, 31), yticks=[10, 25, 40], ylim=(8, 42))
ax.set_xlabel("play", labelpad=1); ax.set_ylabel("mean score, thousands", labelpad=1.5)

S.head_over(axA, "A", "Typing, faster and slower")
S.head_over(axB, "B", "An online game, 30 plays")
out = ROOT / "figures" / "figS4_open_records.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axA], "B": [axB]}, str(ROOT / "figures" / "panels" / "figS4"))
print("slopes", grp, "median intervals", iki.to_dict(), "ratio", round(grp["faster"][0] / grp["slower"][0], 4))
print("game", n_id, share, alpha)
print("written", out)
