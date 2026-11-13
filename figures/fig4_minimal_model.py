"""Figure 4 of the manuscript: the minimal model of practice on a fixed set of 20 elements (code/sim_practice.py).

A, averaged information over 2,000 steps, mean of 300 runs, under random, preferential and fitness-weighted attachment,
which add the same coupling; the right axis gives the time of a 4-bit aiming task from the dynamic law. B, where the
coupling goes under fitness-weighted attachment, in a triangle of shares with the chance hubs (the three strongest of
the other elements at each step) at the top, the peripheral elements (the remaining 14) at the lower left and the fit
elements (0, 1 and 2) at the lower right, drawn as a dotted outline without guide lines, with ticks 0.0, 0.5 and 1.0
on every side. Points every 200 steps give the mean
composition of 300 runs from the even ring of A (shaded marker) and from six other starts in which the peripheral
elements hold the largest share (code/sim_flows.py, sim_out/sim_flows.npz), with an open circle at step 0 and a diamond
at step 2,000. Coloured contour lines give the averaged information, 2.5 to 4 bits, that a composition would have with
coupling even within each group of 3, 14 and 3 elements. The two networks above the triangle are one run from the
even ring at step 0 and at step 2,000, drawn as in Figure 2 and linked by a dotted arrow to the point each shows, the run whose shares at the last step lie closest to the mean (code/sim_networks.py,
sim_out/sim_networks.npz). C, mean degree over 200 runs when a twentieth element is recruited at step 50 with eta = 3,
compared with the targeted and the other elements, with the averaged information of these runs on the right axis.
Palette and type from figures/style.py. Each panel is also written on its own to figures/panels.
Run from the repository root: python figures/fig4_minimal_model.py   (about half a minute)
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(ROOT / "code"))
from style import (Sheet, INK, EDGE, RULE, TONE, GREY, share_colour, tb_axes, save_panels, tern_frame,  # noqa: E402
                   tern_hfield, tern_dots, tern_xy, net_above)
from sim_practice import run  # noqa: E402

N, A_MS, ID = 20, 12.8, 4.0
BETA = 94.7 / np.log2(N)                                             # the first step reproduces Fitts' slope
CEIL = np.log2(N)
RULES = ("random", "pa", "fitness")
H_mean = {rule: np.array([run(rule, seed=s)[0] for s in range(300)]).mean(0) for rule in RULES}
LATE = [run("fitness", steps=3000, seed=s, late=(19, 50, 4.0)) for s in range(200)]
K = np.array([r_[1] for r_ in LATE]); H_late = np.array([r_[0] for r_ in LATE]).mean(0)
top50 = float(np.mean(K[:, -1, 19] > K[:, -1, :19].max(axis=1)))
k_late, k_target, k_other = K[:, :, 19].mean(0), K[:, :, :3].mean((0, 2)), K[:, :, 3:19].mean((0, 2))
total = k_late[-1] + 3 * k_target[-1] + 16 * k_other[-1]
bits = {name: np.log2(k / total * N) for name, k in (("late", k_late[-1]), ("target", k_target[-1]), ("other", k_other[-1]))}
COMP = np.load(ROOT / "sim_out" / "sim_composition.npz")
FL = np.load(ROOT / "sim_out" / "sim_flows.npz")
for rule in ("random", "pa", "fitness"):                             # the composition runs are the runs of panel A
    assert np.max(np.abs(COMP["m4H_" + rule] - H_mean[rule])) < 1e-12

S = Sheet(6.5, 2.38)
axA = S.ax(0.42, 0.40, 1.36, 1.72)
axT = S.ax(2.28, 0.04, 2.00, 2.12)
axC = S.ax(4.74, 0.40, 1.36, 1.72)
NW = np.load(ROOT / "sim_out" / "sim_networks.npz")


# ---------------------------------------------------------------- A, the same coupling, different averaged information
ax = axA
t = np.arange(len(H_mean["fitness"]))
ax.axhline(CEIL, color=RULE, lw=0.6, ls=(0, (2.4, 1.6)))
ax.plot(t, H_mean["random"], color=TONE["small-world"], lw=1.2, ls=(0, (1.0, 1.2)))
ax.plot(t, H_mean["pa"], color=TONE["scale-free"], lw=1.2, ls=(0, (3.0, 1.6)))
ax.plot(t, H_mean["fitness"], color=INK, lw=1.5)
ax.text(1980, CEIL + 0.02, "ceiling log$_2$ 20", ha="right", va="bottom", fontsize=6.0, color=EDGE)
ax.text(1980, H_mean["random"][-1] - 0.03, "random", ha="right", va="top", fontsize=6.2, color=EDGE)
ax.text(1980, H_mean["pa"][-1] - 0.035, "preferential", ha="right", va="top", fontsize=6.2, color=TONE["scale-free"])
ax.text(1980, H_mean["fitness"][-1] - 0.035, "fitness-weighted", ha="right", va="top", fontsize=6.2)
ax.set_xlabel("practice step $t$", labelpad=1); ax.set_ylabel(r"$H(t)$, bits", labelpad=1.5)
ax.set_xlim(0, 2000); ax.set_ylim(3.55, 4.45)
ax2 = ax.twinx(); lo, hi = ax.get_ylim(); ax2.set_ylim(A_MS + BETA * lo * ID, A_MS + BETA * hi * ID)
ax2.set_ylabel("$T$ of a 4-bit task, ms", labelpad=2)
tb_axes(ax, xticks=[0, 1000, 2000], yticks=[3.6, 4.0, 4.4], twin=ax2, twin_ticks=[330, 365, 400])

# ---------------------------------------------------------------- B, where the coupling goes
ax = axT
tern_frame(ax, grid=(), outline="dotted", all_ticks=True, tick_fs=6.4, tick_color=INK, corner_fs=7.2,
           ticks=(0.0, 0.5, 1.0), edge_titles=("chance hubs", "peripheral elements", "fit elements"),
           top=0.30, slide=0.09)                                     # no guide lines, room for the networks above
LAB = [(0.864, 0.112, 0.024), (0.755, 0.202, 0.043), (0.618, 0.315, 0.067), (0.429, 0.471, 0.101)]  # label points
tern_hfield(ax, (3, 14, 3), [2.5, 3.0, 3.5, 4.0], CEIL, fill=False, lw=1.4, label_at=[tuple(np.ravel(tern_xy(c_))) for c_ in LAB])
NS = len(FL["starts"])
for i in range(NS):                                                  # the even ring and six other starts
    tern_dots(ax, FL[f"m4_{i}"][:, [1, 2, 0]], every=200)
x0, y0 = tern_xy(FL["m4_0"][0][[1, 2, 0]])
ax.plot(x0, y0, "o", ms=4.4, color=share_colour(0.0), mec=EDGE, mew=0.6, zorder=7)
ax.text(x0[0], y0[0] - 0.045, "even ring", fontsize=6.4, ha="center", va="top", color=EDGE)
ends = [tern_xy(FL["m4_0"][r_][[1, 2, 0]]) for r_ in (0, -1)]          # where the drawn run starts and ends
ends = [(float(x_[0]), float(y_[0])) for x_, y_ in ends]
NET_W = 0.58                                                         # side of each network, inches
insB = [net_above(S, ax, 0.06, 0.80, NET_W, NW["m4_A0"], "step 0\n$H$ = {:.2f} bits".format(np.log2(20)), "o", ends[0],
                  rad=-0.45, turn=35),
        net_above(S, ax, 0.94, 0.80, NET_W, NW["m4_A"], "step 2,000\n$H$ = {:.2f} bits".format(float(NW["m4_H"])), "D",
                  ends[1], rad=0.20, turn=-75)]

# ---------------------------------------------------------------- C, a late element becomes a hub
ax = axC
t3 = np.arange(len(k_late)); late_plot = k_late.copy(); late_plot[:50] = np.nan
ax.plot(t3, k_other, color=TONE["small-world"], lw=1.2, ls=(0, (1.0, 1.2)))
ax.plot(t3, k_target, color=TONE["scale-free"], lw=1.2, ls=(0, (3.0, 1.6)))
ax.plot(t3, late_plot, color=INK, lw=1.5)
for k_end, name in ((k_late[-1], "late"), (k_target[-1], "target"), (k_other[-1], "other")):
    ax.plot([3000], [k_end], "o", ms=4.0, color=share_colour(bits[name]), mec=EDGE, mew=0.5, zorder=6, clip_on=False)
ax.text(820, 1660, "recruited at step 50,\n" + r"$\eta$ = 3", ha="left", va="top", fontsize=6.2)
axH = ax.twinx()
axH.plot(t3, H_late, color=GREY, lw=1.0, ls=(0, (4.0, 1.4, 1.0, 1.4)), zorder=1)
axH.text(160, H_late[100] + 0.06, "$H(t)$", ha="left", va="bottom", fontsize=6.2, color=GREY)
axH.set_ylabel(r"$H(t)$, bits", labelpad=2)
ax.text(2860, k_target[-1] + 60, r"targeted, $\eta$ = 1.5", ha="right", va="bottom", fontsize=6.2,
        color=TONE["scale-free"])
ax.text(2880, k_other[-1] - 45, r"others, $\eta$ = 1", ha="right", va="top", fontsize=6.2, color=EDGE)
ax.set_xlabel("practice step $t$", labelpad=1); ax.set_ylabel("degree $k$, links", labelpad=1.5)
tb_axes(ax, xticks=[0, 1500, 3000], xlim=(0, 3000), yticks=[0, 900, 1800], ylim=(-60, 1880), twin=axH,
        twin_ticks=[3.0, 3.6, 4.2])
axH.set_ylim(2.75, 4.95)

S.head_over(axA, "A", "Same coupling")
S.head("B", "Where coupling goes", 2.28 + 2.00 / 2)
S.head_over(axC, "C", "A late element becomes a hub")
out = ROOT / "figures" / "fig4_minimal_model.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axA, ax2], "B": [axT] + insB, "C": [axC, axH]}, str(ROOT / "figures" / "panels" / "fig4"))
for rule in H_mean:
    h = H_mean[rule]
    print(f"{rule:8s} H(0) {h[0]:.3f}  H(100) {h[100]:.3f}  H(500) {h[500]:.3f}  H(2000) {h[-1]:.3f}  "
          f"T(0) {A_MS + BETA * h[0] * ID:.1f}  T(2000) {A_MS + BETA * h[-1] * ID:.1f} ms")
for i, st in enumerate(FL["starts"]):
    print(f"panel B, start {tuple(int(x) for x in st)}: shares {np.round(FL[f'm4_{i}'][0], 3)} -> {np.round(FL[f'm4_{i}'][-1], 3)}")
print(f"late recruit at step 50 strongest at step 3000 in {100 * top50:.0f} percent of 200 runs")
print("log2(p N) at step 3000:", {k: round(float(v), 2) for k, v in bits.items()})
print("written", out)
