"""Figure 7 of the manuscript: the Hick–Hyman law as the static case of the dynamic law.

B, the two averaged informations whose product sets the slope (Equation 19): the lights of Hyman's (1953) Experiment II,
Condition 5 (probabilities 1/2, 1/4, 1/8, 1/8, the shares of the hub in Figure 2D), each drawn as a block of width q_j
and height -log2 q_j with its codeword, and even coupling over four elements, as in Figure 2A. A, participant G. C.,
all 24 conditions of Experiments I to III, read from Hyman's Figure 1, marker colour deepening with the stimulus
information (style.demand_colour), with the slope of each experiment, the line Hyman printed for this participant,
RT = 212 + 153 H, and the law in the simulated performer of Supplemental Section S4.4, least-squares lines through its
24 condition means at step 0 (a and beta set from these data) and at step 1,000 (sim_out/sim_static_cases.npz). C, the two kinds of light in Experiment II, Condition 3, as printed by Hyman (p. 194):
the times predicted from the information of each light on his line for G. C., RT = 212 + 153 H, and the times observed.
Palette and type from figures/style.py.
Inputs: data/read_from_figures/hyman1953_fig1_GC_experiment1.csv, data/read_from_figures/hyman1953_fig1_GC_experiments2_3.csv,
data/transcribed/hyman1953_p194_components.csv.
Run from the repository root: python figures/fig7_hyman.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, GREY, EDGE, RULE, TONE, share_colour, outline, demand_colour, tb_axes, save_panels, CLEAR  # noqa: E402

RF, T = ROOT / "data" / "read_from_figures", ROOT / "data" / "transcribed"
e1 = pd.read_csv(RF / "hyman1953_fig1_GC_experiment1.csv")
e23 = pd.read_csv(RF / "hyman1953_fig1_GC_experiments2_3.csv")
comp = pd.read_csv(T / "hyman1953_p194_components.csv").set_index("component")
pts = pd.concat([e1[["experiment", "information_bits", "rt_ms"]], e23[["experiment", "information_bits", "rt_ms"]]])
fit = {e: np.polyfit(g.information_bits, g.rt_ms, 1) for e, g in pts.groupby("experiment")}
bj, aj = np.polyfit(pts.information_bits, pts.rt_ms, 1)
rj = np.corrcoef(pts.information_bits, pts.rt_ms)[0, 1]
assert len(pts) == 24 and abs(bj - 147.8) < 0.2 and abs(rj - 0.985) < 0.001

SIM = np.load(ROOT / "sim_out" / "sim_static_cases.npz")
bl0, al0 = np.polyfit(SIM["x_hyman"], SIM["mean_Hyman_0"].mean(0), 1)
bl1, al1 = np.polyfit(SIM["x_hyman"], SIM["mean_Hyman_1000"].mean(0), 1)
OPEN = "#6E8C99"

S = Sheet(6.5, 1.94)
W_BLK, H_BLK, Y_BLK = 0.70, 1.06, 0.52
XL, XR = 0.38 + 1.92, 1.52 + 1.92                                     # the blocks of panel B
aL = S.ax(XL, Y_BLK, W_BLK, H_BLK)
aR = S.ax(XR, Y_BLK, W_BLK, H_BLK)
axB = S.ax(0.42, 0.42, 1.48, 1.26)                                   # the data, panel A
axC = S.ax(4.84, 0.42, 1.48, 1.26)


def ink_on(rgba):
    r, g, b, _ = rgba
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.45 else INK


def blocks(ax, q, codes, xlabel, hlabel, ylabel=True, label_left=False):
    """Blocks of width q and height -log2 q, drawn in outline in the colour of log2(q n), each with its codeword."""
    q = np.asarray(q, float); n = len(q); I = -np.log2(q); H = float((q * I).sum()); cum = 0.0
    for qi, Ii, cw in zip(q, I, codes):
        c = outline(share_colour(np.log2(qi * n)))
        ax.add_patch(Rectangle((cum, 0), qi, Ii, fc="none", ec=c, lw=1.1, zorder=3 + Ii / 10))
        up = qi < 0.2
        ax.text(cum + qi / 2, min(Ii / 2, (H - 0.15) / 2) if up and Ii > H else Ii / 2, cw, rotation=90 if up else 0,
                fontsize=6.3, ha="center", va="center", color=INK, zorder=7)
        cum += qi
    ax.plot([0, 1], [H, H], color=INK, lw=1.2, zorder=8)
    if label_left:                                                   # clear of the taller blocks on the right
        ax.text(0.02, H + 0.40, hlabel + f" = {H:g}", fontsize=6.3, ha="left", va="bottom", zorder=9,
                bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
    else:
        ax.text(0.97, H + 0.12, hlabel + f" = {H:g}", fontsize=6.4, ha="right", va="bottom", zorder=9,
                bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
    ax.set_xlim(0, 1); ax.set_ylim(0, 4.2); ax.set_yticks([0, 2, 4]); ax.set_xticks([0, 1])
    ax.set_xlabel(xlabel, labelpad=0.8)
    if ylabel:
        ax.set_ylabel("information, bits", labelpad=1.5)
    else:
        ax.set_yticklabels([])
    return H


# ---------------------------------------------------------------- B, the stimulus code and the code of the commands
HS = blocks(aL, [1 / 2, 1 / 4, 1 / 8, 1 / 8], ("0", "10", "110", "111"), r"probabilities $q_j$", "$H_S$",
            label_left=True)
HP = blocks(aR, [1 / 4] * 4, ("00", "01", "10", "11"), r"shares $p_i$", "$H(t)$", ylabel=False)
assert HS == 1.75 and HP == 2.0
PANEL_A_TEXT = [S.text(XL + W_BLK / 2, Y_BLK + H_BLK + 0.07, "lights", fontsize=6.8, ha="center", va="bottom"),
                S.text(XR + W_BLK / 2, Y_BLK + H_BLK + 0.07, "elements", fontsize=6.8, ha="center", va="bottom"),
                S.text((XL + W_BLK + XR) / 2, Y_BLK + H_BLK / 2, r"$\times$", fontsize=11, ha="center", va="center"),
                S.text((XL + XR + W_BLK) / 2, 0.13, r"$RT = a + \beta\,H(t)\,H_S$", fontsize=7.4, ha="center", va="center")]

# ---------------------------------------------------------------- A, three ways of producing the bits, one slope
ax = axB
xx = np.linspace(0, 3.1, 20)
ax.plot(xx, 212 + 153 * xx, color=OPEN, lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
ax.plot(xx, al1 + bl1 * xx, color=TONE["small-world"], lw=0.9, zorder=2)
ax.plot(xx, al0 + bl0 * xx, color=INK, lw=1.0, zorder=2)
STY = {"I": dict(marker="o", ms=3.1), "II": dict(marker="s", ms=2.9), "III": dict(marker="^", ms=3.5)}
LAB = {"I": "number", "II": "probability", "III": "sequence"}
for e, g in pts.groupby("experiment"):
    for xi, yi, ci in zip(g.information_bits, g.rt_ms, demand_colour(g.information_bits, 0, 3)):
        ax.plot(xi, yi, ls="none", zorder=4, mfc=ci, mec=EDGE, mew=0.45, **STY[e])
two = [(1.00, float(e1.loc[e1.alternatives == 2, "rt_ms"].iloc[0])),
       (0.47, float(e23.query("experiment == 'II' and condition == 1").rt_ms.iloc[0]))]
LBL = (0.80, 232)                                          # one label, a leader to each of the two conditions
for xb, yb in two:
    ax.annotate("", xy=(xb, yb), xytext=LBL, arrowprops=dict(arrowstyle="-", lw=0.5, color=RULE, shrinkA=3, shrinkB=2.6))
ax.text(LBL[0], LBL[1], "two lights", fontsize=6.2, ha="center", va="top")
ax.text(-0.05, 748, rf"$RT$ = {aj:.0f} + {bj:.0f} $H_S$", fontsize=6.6, ha="left", **CLEAR)
ax.text(-0.05, 694, f"r = {rj:.3f}".replace("0.", "."), fontsize=6.6, ha="left", style="italic", **CLEAR)
handles = [Line2D([], [], ls="none", mfc="white", mec=EDGE, mew=0.6, **STY[e]) for e in ("I", "II", "III")]
leg = ax.legend(handles, [f"{LAB[e]}, {fit[e][0]:.0f}" for e in ("I", "II", "III")], loc="upper left",
                bbox_to_anchor=(-0.02, 0.80), fontsize=5.9, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.12,
                handlelength=1.0, title="slope, ms per bit", title_fontsize=5.9)
ax.add_artist(leg)
ax.legend([Line2D([], [], color=OPEN, lw=0.9, ls=(0, (3.0, 1.6))), Line2D([], [], color=INK, lw=1.0),
           Line2D([], [], color=TONE["small-world"], lw=0.9)], ["Hyman’s line", "law, step 0", "law, step 1,000"],
          loc="lower right", fontsize=5.9, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.12, handlelength=1.4)
tb_axes(ax, xticks=[0, 1.5, 3], xlim=(-0.15, 3.2), yticks=[200, 450, 700], ylim=(150, 800))
ax.set_xlabel(r"stimulus information $H_S$, bits", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)

# ---------------------------------------------------------------- C, each light and the mean of the condition
ax = axC
xx = np.linspace(0, 4.4, 20)
ax.plot(xx, 212 + 153 * xx, color=OPEN, lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
rows = [("frequent light", 13 / 16, "13/16"), ("each rare light", 3 / 16, "3/16 in all")]
for k, w, lab in rows:
    x_, yo, yp = comp.loc[k, ["information_bits", "rt_observed_ms", "rt_predicted_ms"]]
    c = share_colour(np.log2((w if w > 0.5 else w / 3) * 4))
    ax.plot(x_, yp, "o", ms=3.6, color="white", mec=INK, mew=0.7, zorder=4)
    ax.annotate("", xy=(x_, yo), xytext=(x_, yp), arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.7,
                color=EDGE, shrinkA=3, shrinkB=4), zorder=3)
    ax.plot(x_, yo, "o", ms=3.4 + 2.4 * w, color=c, mec=EDGE, mew=0.5, zorder=5)
    if w > 0.5:
        ax.text(x_, yo + 34, lab, fontsize=6.3, ha="center", va="bottom")
    else:
        ax.text(x_ - 0.22, yo, lab, fontsize=6.3, ha="right", va="center")
xm, ym, pm = comp.loc["condition mean", ["information_bits", "rt_observed_ms", "rt_predicted_ms"]]
ax.plot(xm, pm, "D", ms=3.6, color="white", mec=INK, mew=0.7, zorder=4)
ax.plot(xm, ym, "D", ms=2.4, color=INK, zorder=6)
ax.annotate("condition mean", xy=(xm + 0.05, ym - 8), xytext=(1.35, 300), fontsize=6.3, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=RULE, shrinkA=1, shrinkB=1))
ax.text(-0.05, 932, r"$RT$ = 212 + 153 $H$", fontsize=6.6, ha="left", **CLEAR)
ax.legend([Line2D([], [], ls="none", marker="o", ms=3.6, color="white", mec=INK, mew=0.7),
           Line2D([], [], ls="none", marker="o", ms=3.6, color=GREY, mec=EDGE, mew=0.5)],
          ["from its information", "observed"], loc="upper left", bbox_to_anchor=(-0.02, 0.89), fontsize=6.2, handletextpad=0.2, borderaxespad=0.1,
          labelspacing=0.15, handlelength=1.0)
tb_axes(ax, xticks=[0, 2, 4], xlim=(-0.15, 4.4), yticks=[200, 550, 900], ylim=(150, 1000))
ax.set_xlabel(r"information of a light, $-\log_2 q_j$, bits", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)

S.head("B", "Two averaged informations", (XL + XR + W_BLK) / 2)
S.head_over(axB, "A", "Three ways, one slope")
S.head_over(axC, "C", "Each light and the mean")
out = ROOT / "figures" / "fig7_hyman.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axB], "B": [aL, aR] + PANEL_A_TEXT, "C": [axC]}, str(ROOT / "figures" / "panels" / "fig7"))
print(f"law in the simulated performer: step 0 RT = {al0:.1f} + {bl0:.1f} H, step 1,000 RT = {al1:.1f} + {bl1:.1f} H")
print(f"all 24: RT = {aj:.1f} + {bj:.1f} H, r = {rj:.4f}; slopes", {e: round(float(v[0]), 1) for e, v in fit.items()},
      "intercepts", {e: round(float(v[1]), 1) for e, v in fit.items()})
print("two lights:", two)
print("written", out)
