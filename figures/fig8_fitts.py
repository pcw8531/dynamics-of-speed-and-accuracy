"""Figure 8 of the manuscript: Fitts' law and Woodworth's two components as the static case of the dynamic law.

B, a scheme of one aimed movement (Equations 20 and 21). The initial adjustment brings the remaining error from
log2(2D/W) bits to the level it leaves, and current control then removes log2(1/rho) bits with each correction of
duration tau = beta L log2(1/rho), until the hand is within W/2 of the centre of the target. With the eyes shut no
correction is made and the error stays where the initial adjustment left it. After practice, H(t) is lower, so each
correction is shorter and the movement ends sooner, by more the more corrections it needs. Drawn to scale for six
bits, one bit per correction and H falling from 4.322 to 3.005 bits as in Figure 2. A, Fitts (1954) Table 1, mean
movement time of 16 participants for the one-ounce and one-pound styluses, one-ounce markers deepening in colour with
ID (style.demand_colour), with Fitts' own reading of his table, a constant rate MT = ID / I_P at the mean of the 16
printed indices of performance for the one-ounce stylus, the least-squares line of MacKenzie (1992), and the law in the
simulated performer of Supplemental Section S4.4 at step 0 (a0 = 180 ms and beta set from this table) and at step 1,000
(sim_out/sim_static_cases.npz). C, the pairs (log2(1/rho), tau) that give one slope b = tau / log2(1/rho), for the slope of A and for the slope
after the same practice. The point at two bits per correction is an illustration. Palette and type from figures/style.py.
Input: data/transcribed/fitts1954_table1.csv.
Run from the repository root: python figures/fig8_fitts.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, GREY, EDGE, RULE, TRACK, TONE, demand_colour, tb_axes, save_panels, CLEAR  # noqa: E402
sys.path.insert(0, str(ROOT / "code"))
import sim_full as sf  # noqa: E402

f = pd.read_csv(ROOT / "data" / "transcribed" / "fitts1954_table1.csv")
b, a = np.polyfit(f.Id_bits, 1000 * f.t_1oz_s, 1)
res = 1000 * f.t_1oz_s - (a + b * f.Id_bits)
H0 = np.log2(20)                                             # even coupling over 20 elements
HT = float(np.mean([sf.run(N=20, n_start=20, steps=1000, lam=0.5, t_elim=200, seed=s)[0][-1] for s in range(300)]))
assert abs(HT - 3.0055) < 1e-3                               # Figure 2B at step 1,000, mean of 300 runs
B_LAW = 94.7                                                 # the slope as used in Figure 2 and the simulations
bT = B_LAW * HT / H0
assert abs(a - 12.8) < 0.05 and abs(b - 94.7) < 0.05

S = Sheet(6.5, 1.94)
axA = S.ax(2.42, 0.42, 1.70, 1.26)                                   # the scheme, panel B
axB = S.ax(0.42, 0.42, 1.48, 1.26)                                   # the data, panel A
axC = S.ax(4.84, 0.42, 1.48, 1.26)
BEG, PRA = TONE["small-world"], INK

# ---------------------------------------------------------------- B, initial adjustment and current control
ax = axA
ID, Y0, T0, TAU = 6.0, 3.0, 1.15, 1.0
TAU_T = TAU * HT / H0
ax.axhspan(-0.6, 0.0, color=TRACK, lw=0, zorder=0)
ax.text(0.08, -0.3, "within $W/2$ of the centre", fontsize=6.0, ha="left", va="center", color=GREY)
s_ = np.linspace(0, 1, 40)
ax.plot(T0 * s_, Y0 + (ID - Y0) * (1 + np.cos(np.pi * s_)) / 2, color=INK, lw=1.2, zorder=4)


def stairs(tau, n_steps):
    xs, ys, t, y = [T0], [Y0], T0, Y0
    for _ in range(n_steps):
        t += tau; xs += [t, t]; ys += [y, y - 1.0]; y -= 1.0
    return np.array(xs), np.array(ys)


xb, yb = stairs(TAU, 3)
xp, yp = stairs(TAU_T, 3)
ax.plot([T0, 5.0], [Y0, Y0], color=TONE["scale-free"], lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
ax.plot(xb, yb, color=BEG, lw=1.3, zorder=3, solid_joinstyle="miter")
ax.plot(xp, yp, color=PRA, lw=1.1, zorder=4, solid_joinstyle="miter")
ax.text(4.95, Y0 + 0.18, "eyes shut", fontsize=6.2, ha="right", va="bottom", color=TONE["scale-free"])
for x_, lab, c in ((xb[-1], "before", BEG), (xp[-1], "after practice", PRA)):
    ax.plot([x_, x_], [-0.6, 0.0], color=c, lw=0.8)
ax.text(xb[-1] + 0.06, 0.55, "before", fontsize=6.2, ha="left", va="center", color=TONE["scale-free"])
ax.text(xp[-1] - 0.08, 0.36, "after\npractice", fontsize=6.0, ha="right", va="center", linespacing=0.92,
        bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"), zorder=5)
ax.annotate("", xy=(T0, 3.55), xytext=(T0 + TAU, 3.55), arrowprops=dict(arrowstyle="<->", lw=0.5, color=EDGE,
            mutation_scale=4, shrinkA=0, shrinkB=0))
ax.text(T0 + TAU / 2, 3.68, r"$\tau$", fontsize=7.0, ha="center", va="bottom")
ax.annotate("", xy=(T0 + TAU + 0.12, Y0), xytext=(T0 + TAU + 0.12, Y0 - 1.0), arrowprops=dict(arrowstyle="<->",
            lw=0.5, color=EDGE, mutation_scale=4, shrinkA=0, shrinkB=0))
ax.text(T0 + TAU + 0.2, Y0 - 0.5, r"$\log_2(1/\rho)$", fontsize=6.6, ha="left", va="center")
ax.text(T0 / 2, ID + 0.38, "initial\nadjustment", fontsize=6.2, ha="left", va="bottom")
ax.text(T0 + 1.75, 4.15, "current control", fontsize=6.2, ha="left", va="center")
ax.set_xlim(0, 5.0); ax.set_ylim(-0.6, 7.6); ax.set_yticks([0, 6]); ax.set_yticklabels(["0", r"$ID$"])
ax.set_xticks([]); ax.set_xlabel("time", labelpad=2); tb_axes(ax, grid=False)
ax.set_ylabel(r"error left, $\log_2(2E/W)$, bits", labelpad=1.5)

# ---------------------------------------------------------------- A, Fitts (1954), Table 1
ax = axB
SIM = np.load(ROOT / "sim_out" / "sim_static_cases.npz")
xs, law0, law1 = SIM["x_fitts"], SIM["mean_Fitts_0"].mean(0), SIM["mean_Fitts_1000"].mean(0)
IP = float(f.Ip_1oz_bits_per_s.mean())                              # mean of the 16 printed indices of performance
OPEN = "#6E8C99"
xx = np.linspace(0.6, 7.4, 20)
ax.plot(xx, 1000 * xx / IP, color=GREY, lw=0.8, ls=(0, (1.0, 1.3)), zorder=2)
ax.plot(xx, a + b * xx, color=OPEN, lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
ax.plot(xs, law1, color=TONE["small-world"], lw=0.9, zorder=2)
ax.plot(xs, law0, color=INK, lw=1.0, zorder=2)
ax.plot(f.Id_bits, 1000 * f.t_1lb_s, "s", ms=2.6, mfc="white", mec=OPEN, mew=0.7, zorder=3)
for xi, yi, ci in zip(f.Id_bits, 1000 * f.t_1oz_s, demand_colour(f.Id_bits, 1, 7)):
    ax.plot(xi, yi, "o", ms=3.0, mfc=ci, mec=EDGE, mew=0.45, zorder=4)
y_id1 = float(1000 * f.loc[f.Id_bits == 1, "t_1oz_s"].iloc[0])
ax.annotate("$W = D$,\nno correction", xy=(1.0, y_id1 + 14), xytext=(0.62, 330), fontsize=6.0, ha="left", va="bottom", **CLEAR,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=RULE, shrinkA=1, shrinkB=2))
ax.text(0.7, 815, rf"$MT$ = {a:.1f} + {b:.1f} $ID$", fontsize=6.6, ha="left", **CLEAR)
leg = ax.legend([Line2D([], [], ls="", marker="o", ms=3.0, mfc=demand_colour(4, 1, 7), mec=EDGE, mew=0.45),
                 Line2D([], [], ls="", marker="s", ms=2.6, mfc="white", mec=OPEN, mew=0.7)],
                ["1-oz stylus", "1-lb stylus"], loc="upper left", bbox_to_anchor=(-0.02, 0.90), fontsize=5.9,
                handletextpad=0.2, borderaxespad=0.1, labelspacing=0.12, handlelength=1.0)
ax.add_artist(leg)
ax.legend([Line2D([], [], color=GREY, lw=0.8, ls=(0, (1.0, 1.3))), Line2D([], [], color=OPEN, lw=0.9, ls=(0, (3.0, 1.6))),
           Line2D([], [], color=INK, lw=1.0), Line2D([], [], color=TONE["small-world"], lw=0.9)],
          [f"Fitts’ rate, {IP:.1f} bits/s", "MacKenzie’s line", "law, step 0", "law, step 1,000"], loc="lower right",
          fontsize=5.9, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.12, handlelength=1.4)
ax.set_xlabel(r"$ID = \log_2(2D/W)$, bits", labelpad=1); ax.set_ylabel("movement time, ms", labelpad=1.5)
tb_axes(ax, xticks=[1, 4, 7], xlim=(0.5, 7.5), yticks=[0, 400, 800], ylim=(0, 880))

# ---------------------------------------------------------------- C, two ways to a lower slope
ax = axC
xx = np.linspace(0, 3.2, 20)
ax.plot(xx, B_LAW * xx, color=BEG, lw=1.3)
ax.plot(xx, bT * xx, color=PRA, lw=1.1)
P = (2.0, B_LAW * 2.0)
ax.annotate("", xy=(2.0, bT * 2.0), xytext=P, arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.8, color=EDGE,
            shrinkA=3, shrinkB=2))
ax.annotate("", xy=(B_LAW * 2.0 / bT, B_LAW * 2.0), xytext=P, arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.8,
            color=EDGE, shrinkA=3, shrinkB=2))
ax.plot(*P, "o", ms=3.4, color="white", mec=EDGE, mew=0.7, zorder=5)
ax.annotate("shorter\ncorrection", xy=(2.02, (B_LAW + bT) * 1.0), xytext=(2.30, 52), fontsize=6.2, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=RULE, shrinkA=1, shrinkB=1))
ax.text((2.0 + B_LAW * 2.0 / bT) / 2, B_LAW * 2.0 + 9, "larger\ncorrection", fontsize=6.2, ha="center", va="bottom")
ax.text(0.95, B_LAW * 0.95 + 24, f"$b$ = {B_LAW:.1f}", fontsize=6.6, ha="right", va="bottom", color=TONE["scale-free"])
ax.text(2.62, bT * 2.62 - 22, f"$b$ = {bT:.1f}", fontsize=6.6, ha="left", va="top")
tb_axes(ax, xticks=[0, 1.5, 3], xlim=(0, 3.2), yticks=[0, 150, 300], ylim=(0, 310))
ax.set_xlabel(r"bits per correction, $\log_2(1/\rho)$", labelpad=1); ax.set_ylabel(r"$\tau$, ms per correction", labelpad=1.5)

S.head_over(axA, "B", "Woodworth’s two components")
S.head_over(axB, "A", "Fitts (1954)")
S.head_over(axC, "C", "Two ways to a lower slope")
out = ROOT / "figures" / "fig8_fitts.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axB], "B": [axA], "C": [axC]}, str(ROOT / "figures" / "panels" / "fig8"))
print(f"1-oz: MT = {a:.1f} + {b:.1f} ID; residual by ID", res.groupby(f.Id_bits).mean().round(1).to_dict())
print(f"slope after practice {bT:.2f} ms per bit; at 2 bits per correction tau {B_LAW * 2:.1f} -> {bT * 2:.1f} ms, or "
      f"{B_LAW * 2 / bT:.2f} bits per correction at tau {B_LAW * 2:.1f} ms")
print("written", out)
