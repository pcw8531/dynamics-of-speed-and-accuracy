"""Figure 3 of the manuscript: one practice step of the simulations as a flow chart, read from left to right.

The network of N elements is built once, with the links A_ij between elements, the degree k_i and the fitness
eta_i = 1 + lambda g_i of each element and H(0) = log2 N for the even ring. Each step then recruits an element when one
is due, and otherwise forms one link and, in the full model from step 200, eliminates one, before the shares and H(t)
are updated and the law gives the slope and the time of the act. The number beside each box gives the lines of
Table 1 of the manuscript that carry out its operation, and Section S5.1 of the Supplemental Material lists the code of
the step, condensed from the function run of code/sim_full.py. Drawn in black only. No data enter this figure.
Run from the repository root: python figures/fig3_practice_step.py
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib.patches import FancyBboxPatch, Polygon, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet  # noqa: E402

K = "#000000"                                                      # every line and letter of the chart in black
FW, FH = 6.5, 1.62
S = Sheet(FW, FH)
ax = S.ax(0.0, 0.0, FW, FH); ax.set_xlim(0, FW); ax.set_ylim(0, FH); ax.axis("off")
FS, FK = 6.2, 5.6                                                  # text of the boxes, line numbers and branch labels


# ---------------------------------------------------------------- shapes
def box(xc, yc, w, h, lines, key=None):
    ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, fc="white", ec=K, lw=0.6,
                                boxstyle="round,pad=0,rounding_size=0.04", zorder=3))
    n = len(lines)
    for j, ln in enumerate(lines):
        ax.text(xc, yc + (n - 1) * 0.052 - j * 0.104, ln, fontsize=FS, ha="center", va="center", color=K, zorder=4)
    if key:
        ax.text(xc - w / 2 + 0.01, yc + h / 2 + 0.025, key, fontsize=FK, ha="left", va="bottom", color=K, zorder=4)


def diamond(xc, yc, w, h, lines, key=None, key_side="left"):
    ax.add_patch(Polygon([(xc - w / 2, yc), (xc, yc + h / 2), (xc + w / 2, yc), (xc, yc - h / 2)], closed=True,
                         fc="white", ec=K, lw=0.6, zorder=3))
    n = len(lines)
    for j, ln in enumerate(lines):
        ax.text(xc, yc + (n - 1) * 0.048 - j * 0.096, ln, fontsize=FS - 0.2, ha="center", va="center", color=K, zorder=4)
    if key:
        dx = -0.05 if key_side == "left" else 0.05
        ax.text(xc + dx, yc + h / 2 + 0.005, key, fontsize=FK, ha="right" if key_side == "left" else "left",
                va="bottom", color=K, zorder=4)


def arrow(pts, label=None, at=None, ha="center"):
    """A path of straight segments with an arrowhead at the last point."""
    pts = [np.asarray(p, float) for p in pts]
    for p0, p1 in zip(pts[:-2], pts[1:-1]):
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=K, lw=0.6, zorder=2, solid_capstyle="butt")
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=5.5, lw=0.6, color=K,
                                 shrinkA=0, shrinkB=0, zorder=2))
    if label:
        ax.text(*at, label, fontsize=FK, ha=ha, va="center", style="italic", color=K, zorder=4)


# ---------------------------------------------------------------- the chart
M, LOOP, SIDE, LOW = 0.98, 1.50, 0.38, 0.08                         # main row, loop line, side row, lowest return
X = dict(net=0.56, step=1.36, rec=2.04, form=2.85, elim=3.67, share=4.45, law=5.26, last=6.07)
W = dict(net=0.92, step=0.46, rec=0.64, form=0.72, elim=0.66, share=0.56, law=0.80, last=0.62, side=0.74, out=0.58)
H_BOX, H_DIA, H_NET = 0.32, 0.42, 0.56

box(X["net"], M, W["net"], H_NET, ["network at step 0", r"$N$ elements, links $A_{ij}$", r"$k_i$ and $\eta_i$",
                                    r"$H(0) = \log_2 N$"], key="1 to 3")
box(X["step"], M, W["step"], 0.30, ["practice", r"step $t$"], key="4")
diamond(X["rec"], M, W["rec"], H_DIA, ["recruit", r"at $t$?"], key="5")
box(X["rec"], SIDE, W["side"], H_BOX, ["recruitment", r"$m$ links, by $\eta_j k_j$"], key="6")
box(X["form"], M, W["form"], H_BOX, ["formation", r"one link, by $\Pi_i$"], key="7, 8")
diamond(X["elim"], M, W["elim"], H_DIA, [r"$t \geq 200$,", "full model?"], key="9")
box(X["elim"], SIDE, W["side"], H_BOX, ["elimination", r"one end by $1/\eta_i$"], key="10")
box(X["share"], M, W["share"], H_BOX, [r"shares $p_i$", r"and $H(t)$"], key="11")
box(X["law"], M, W["law"], H_BOX, [r"$b(t) = \beta H(t)$", r"$T = a + \beta H(t)\,\log_2 X$"], key="12")
diamond(X["last"], M, W["last"], H_DIA, [r"$t$ < last", "step?"], key="4")
box(X["last"], SIDE, W["out"], H_BOX, ["output", r"$H(t)$ and $T(t)$"], key="13")

e = lambda name, w=None: (X[name] - (w or W[name]) / 2, X[name] + (w or W[name]) / 2)     # left and right edges
arrow([(e("net")[1], M), (e("step")[0], M)])
arrow([(e("step")[1], M), (e("rec")[0], M)])
arrow([(e("rec")[1], M), (e("form")[0], M)], "no", ((e("rec")[1] + e("form")[0]) / 2, M + 0.06))
arrow([(X["rec"], M - H_DIA / 2), (X["rec"], SIDE + H_BOX / 2)], "yes", (X["rec"] + 0.04, (M - H_DIA / 2 + SIDE + H_BOX / 2) / 2), "left")
arrow([(e("form")[1], M), (e("elim")[0], M)])
arrow([(e("elim")[1], M), (e("share")[0], M)], "no", ((e("elim")[1] + e("share")[0]) / 2, M + 0.06))
arrow([(X["elim"], M - H_DIA / 2), (X["elim"], SIDE + H_BOX / 2)], "yes", (X["elim"] + 0.04, (M - H_DIA / 2 + SIDE + H_BOX / 2) / 2), "left")
XE, XR = X["share"] - 0.12, X["share"] + 0.12                       # where the two side paths enter the shares
arrow([(e("elim", W["side"])[1], SIDE), (XE, SIDE), (XE, M - H_BOX / 2)])
arrow([(X["rec"], SIDE - H_BOX / 2), (X["rec"], LOW), (XR, LOW), (XR, M - H_BOX / 2)])
arrow([(e("share")[1], M), (e("law")[0], M)])
arrow([(e("law")[1], M), (e("last")[0], M)])
arrow([(X["last"], M + H_DIA / 2), (X["last"], LOOP), (X["step"], LOOP), (X["step"], M + 0.15)], "yes",
      (X["last"] + 0.04, (M + H_DIA / 2 + LOOP) / 2), "left")
arrow([(X["last"], M - H_DIA / 2), (X["last"], SIDE + H_BOX / 2)], "no", (X["last"] + 0.04, (M - H_DIA / 2 + SIDE + H_BOX / 2) / 2), "left")

out = ROOT / "figures" / "fig3_practice_step.png"
S.fig.savefig(out)
print("written", out)
