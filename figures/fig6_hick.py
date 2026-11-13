"""Figure 6 of the manuscript: Hick's law as the static case of the dynamic law.

B, a choice among n = 8 equally likely lights made by three halvings, each a command coded over the performer's
elements, so RT = a + beta (L_1 + L_2 + L_3) >= a + 3 beta H(t) (Equation 18). A, participant G. C. of Hyman (1953),
Experiment I, one to eight lights, read from his Figure 1, marker colour deepening with the bits of the choice
(style.demand_colour), with the line Hyman printed for this participant, RT = 212 + 153 H, and the law in the simulated
performer of Supplemental Section S4.4 at step 0 (a and beta set from these data) and at step 1,000
(sim_out/sim_static_cases.npz). C, the pairs (H, beta) that
give the slope of A: even coupling over 4 or 20 elements, and practice at a fixed beta, which lowers H from 4.322 to
3.005 bits as in Figure 2 and the slope by the same ratio (Equation 10). Palette and type from figures/style.py.
Input: data/read_from_figures/hyman1953_fig1_GC_experiment1.csv.
Run from the repository root: python figures/fig6_hick.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, GREY, EDGE, LINK, RULE, TONE, share_colour, demand_colour, tb_axes, save_panels, CLEAR  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
sys.path.insert(0, str(ROOT / "code"))
import sim_full as sf  # noqa: E402

e1 = pd.read_csv(ROOT / "data" / "read_from_figures" / "hyman1953_fig1_GC_experiment1.csv")
x, y = np.log2(e1.alternatives.values), e1.rt_ms.values
b, a = np.polyfit(x, y, 1)
r = np.corrcoef(x, y)[0, 1]
H20, H4 = np.log2(20), 2.0                                 # even coupling over 20 and over 4 elements
H_T = float(np.mean([sf.run(N=20, n_start=20, steps=1000, lam=0.5, t_elim=200, seed=s)[0][-1] for s in range(300)]))
assert abs(H_T - 3.0055) < 1e-3                             # Figure 2B at step 1,000, mean of 300 runs
assert abs(b - 147.8) < 0.1 and abs(a - 226.6) < 0.2 and abs(r - 0.993) < 0.001

SIM = np.load(ROOT / "sim_out" / "sim_static_cases.npz")
xs = SIM["x_hick"]; law0 = SIM["mean_Hick_0"].mean(0); law1 = SIM["mean_Hick_1000"].mean(0)
assert np.allclose(SIM["obs_hick"], y, atol=0.05)                   # the simulated performer used these same values
OPEN = "#6E8C99"

S = Sheet(6.5, 2.0)
axA = S.ax(2.04, 0.16, 2.30, 1.56)                                   # the scheme, panel B
axB = S.ax(0.42, 0.40, 1.44, 1.22)                                   # the data, panel A
axC = S.ax(4.86, 0.40, 1.44, 1.22)

# ---------------------------------------------------------------- B, halvings and commands
ax = axA
ax.set_xlim(0, 10); ax.set_ylim(0, 6.8); ax.axis("off")
LEV = [5.95, 4.55, 3.15, 1.75]                              # root, three halvings down to the eight lights
XL = np.linspace(0.45, 5.55, 8)
nodes = [[XL.reshape(2 ** d, -1).mean(1)] for d in range(4)]
pick = 5                                                    # the light shown, 101 in binary
path = [pick >> (3 - d) for d in range(4)]
for d in range(3):
    for j, xp in enumerate(nodes[d][0]):
        for c in (0, 1):
            k = 2 * j + c; xc = nodes[d + 1][0][k]
            on = path[d] == j and path[d + 1] == k
            ax.plot([xp, xc], [LEV[d], LEV[d + 1]], color=INK if on else LINK, lw=1.2 if on else 0.6,
                    zorder=3 if on else 2, solid_capstyle="round")
            if on:
                ax.text((xp + xc) / 2 + (0.16 if c else -0.16), (LEV[d] + LEV[d + 1]) / 2, str(c), fontsize=6.4,
                        ha="left" if c else "right", va="center")
for d in range(3):
    for xp in nodes[d][0]:
        ax.add_patch(Circle((xp, LEV[d]), 0.09, fc=EDGE, ec="none", zorder=4))
for k, xp in enumerate(XL):
    ax.add_patch(Circle((xp, LEV[3]), 0.24, fc=INK if k == pick else "white", ec=EDGE, lw=0.6, zorder=4))
for d in range(3):                                          # one command for each halving
    ym = (LEV[d] + LEV[d + 1]) / 2
    ax.text(6.05, ym, f"halving {d + 1}", fontsize=6.6, ha="left", va="center")
    ax.text(7.85, ym, rf"$\beta L_{d + 1}$", fontsize=7.2, ha="left", va="center")
ax.plot([7.62, 7.62], [LEV[3] + 0.32, LEV[0] - 0.05], color=RULE, lw=0.6)
ax.text(5.0, 0.72, r"$RT = a + \beta\,(L_1 + L_2 + L_3) \geq a + \beta H(t)\,\log_2 8$", fontsize=7.2, ha="center",
        va="center")
ax.text(7.85, LEV[0] + 0.02, "time", fontsize=6.4, ha="left", va="center", color=GREY)

# ---------------------------------------------------------------- A, Hyman (1953), Experiment I
ax = axB
xx = np.linspace(0, 3.1, 20)
ax.plot(xx, 212 + 153 * xx, color=OPEN, lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
ax.plot(xs, law1, color=TONE["small-world"], lw=0.9, zorder=2)
ax.plot(xs, law0, color=INK, lw=1.0, zorder=2)
for xi, yi, ci in zip(x, y, demand_colour(x, 0, 3)):
    ax.plot(xi, yi, "o", ms=3.2, mfc=ci, mec=EDGE, mew=0.45, zorder=3)
ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=3.0, mfc=demand_colour(1.5, 0, 3), mec=EDGE, mew=0.45),
                   Line2D([], [], color=OPEN, lw=0.9, ls=(0, (3.0, 1.6))), Line2D([], [], color=INK, lw=1.0),
                   Line2D([], [], color=TONE["small-world"], lw=0.9)],
          labels=["G. C., observed", "Hyman’s line", "law, step 0", "law, step 1,000"], loc="upper left",
          bbox_to_anchor=(-0.02, 0.83),
          fontsize=5.9, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.12, handlelength=1.4)
ax.text(-0.05, 742, rf"$RT$ = {a:.0f} + {b:.0f} $\log_2 n$", fontsize=6.6, ha="left", **CLEAR)
ax.text(-0.05, 688, f"r = {r:.3f}".replace("0.", "."), fontsize=6.6, ha="left", style="italic", **CLEAR)
tb_axes(ax, xticks=[0, 1.5, 3], xlim=(-0.15, 3.2), yticks=[200, 450, 700], ylim=(150, 800))
ax.set_xlabel(r"$\log_2 n$, bits", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)
top = ax.secondary_xaxis("top"); top.set_xticks(np.log2([1, 2, 4, 8])); top.set_xticklabels(["1", "2", "4", "8"])
top.tick_params(length=1.6, pad=1, labelsize=6.4); top.spines["top"].set_linewidth(0.6)
ax.text(3.32, 800, "$n$", fontsize=6.8, ha="left", va="center", clip_on=False)

# ---------------------------------------------------------------- C, one slope fixes only the product
ax = axC
hh = np.linspace(1.0, 5.0, 200)
b_t = b * H_T / H20
ax.plot(hh, b / hh, color=INK, lw=1.1)
ax.plot(hh, b_t / hh, color=TONE["small-world"], lw=1.1, ls=(0, (3.0, 1.6)))
for H_, n_, dx, ha in ((H4, 4, 0.1, "left"), (H20, 20, 0.05, "center")):
    ax.plot(H_, b / H_, "o", ms=4.0, color=share_colour(0.0), mec=EDGE, mew=0.6, zorder=4)
    ax.text(H_ + dx, b / H_ + 7, rf"$N$ = {n_}, $\beta$ = {b / H_:.0f}", fontsize=6.4, ha=ha, va="bottom")
ax.annotate("", xy=(H_T + 0.05, b / H20), xytext=(H20 - 0.09, b / H20),
            arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.8, color=INK), zorder=5)
ax.plot(H_T, b / H20, "o", ms=3.4, color=INK, mec="white", mew=0.4, zorder=6)
ax.text((H_T + H20) / 2, b / H20 - 9, "practice", fontsize=6.4, ha="center", va="top")
ax.text(1.12, 141, f"$b$ = {b:.0f}", fontsize=6.6, ha="left", va="center")
ax.text(1.1, 46, f"$b$ = {b_t:.0f}", fontsize=6.6, ha="left", va="center", color=TONE["scale-free"])
tb_axes(ax, xticks=[1, 3, 5], xlim=(1.0, 5.0), yticks=[0, 75, 150], ylim=(0, 160))
ax.set_xlabel(r"averaged information $H(t)$, bits", labelpad=1); ax.set_ylabel(r"$\beta$, ms per command bit", labelpad=1.5)

S.head("B", "Halvings and commands", 2.04 + 2.30 / 2)
S.head_over(axB, "A", "Hyman (1953), Experiment I")
S.head_over(axC, "C", r"One slope, many $\beta$ and $H$")
out = ROOT / "figures" / "fig6_hick.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axB, top], "B": [axA], "C": [axC]}, str(ROOT / "figures" / "panels" / "fig6"))
print(f"Experiment I: RT = {a:.1f} + {b:.1f} log2 n, r = {r:.4f}; beta = {b / H4:.1f} (N = 4), {b / H20:.2f} (N = 20); "
      f"slope after practice {b_t:.1f} ms per bit")
print("written", out)
