"""Figure S2 of the Supplemental Material: the classical laws as static cases, by simulation (Section S4.4).

One simulated performer, coupling of the full model at step 0 (even over 20 elements) and at step 1,000 (Figure 5),
commands coded three at a time over the shares of coupling (code/sim_static_cases.py). A, Hick: choices among 1 to 8
equally likely lights, with the times of participant G. C. (Hyman, 1953, Experiment I). B, Hyman: the 24 conditions of
G. C., condition means. C, Fitts: an initial adjustment and corrections of one bit, with Fitts (1954, Table 1, 1-oz
stylus). Simulated points are means of 200 repetitions of 3,000 trials; lines are least-squares fits (C, ID 2 to 7).
D, the slope at step 1,000 divided by the slope at step 0 in each task, mean and 95 percent of repetitions, against the
ratio of averaged information. Palette and type from figures/style.py.
Run from the repository root after code/sim_static_cases.py: python figures/figS2_static_cases.py sim_out
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, GREY, EDGE, RULE, TONE, CLEAR, tb_axes, save_panels  # noqa: E402

d = np.load(Path(sys.argv[1] if len(sys.argv) > 1 else "sim_out") / "sim_static_cases.npz")
H0, H1, L0, L1 = (float(d[k]) for k in ("H0", "H1000", "L0", "L1000"))
BEG, PRA, OBS = TONE["small-world"], INK, TONE["scale-free"]

S = Sheet(6.5, 1.94)
axes = [S.ax(x, 0.42, 1.10, 1.26) for x in (0.42, 1.98, 3.54)]
axD = S.ax(5.22, 0.42, 1.10, 1.26)


def task(ax, x, key, use=None, obs=None, x_obs=None):
    use = np.ones(len(x), bool) if use is None else use
    xx = np.linspace(x[use].min(), x[use].max(), 20)
    out = {}
    for st, c, z in ((0, BEG, 3), (1000, PRA, 4)):
        m = d[f"mean_{key}_{st}"].mean(0)
        b, a = np.polyfit(x[use], m[use], 1); out[st] = (a, b)
        ax.plot(xx, a + b * xx, color=c, lw=0.9, zorder=2)
        ax.plot(x, m, "o", ms=2.7, color=c, mec="white", mew=0.3, zorder=z)
    if obs is not None:
        ax.plot(x_obs, obs, "s", ms=2.9, mfc="none", mec=OBS, mew=0.7, zorder=5)
    return out


def slope_labels(ax, f, x0, x1, up, down, fmt):
    """Slope of step 0 above its line, ending at x0; slope of step 1,000 below its line, starting at x1."""
    ax.text(x0, f[0][0] + f[0][1] * x0 + up, fmt.format(f[0][1]) + " ms per bit", fontsize=6.3, ha="right",
            va="bottom", color=TONE["random"])
    ax.text(x1, f[1000][0] + f[1000][1] * x1 - down, fmt.format(f[1000][1]) + " ms per bit", fontsize=6.3,
            ha="left", va="top")


ax = axes[0]
f = task(ax, d["x_hick"], "Hick", obs=d["obs_hick"], x_obs=d["x_hick"])
tb_axes(ax, xticks=[0, 1.5, 3], xlim=(-0.15, 3.15), yticks=[200, 450, 700], ylim=(150, 780))
ax.set_xlabel(r"$\log_2 n$, bits", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)
slope_labels(ax, f, 2.05, 1.45, 30, 75, "{:.0f}")
ax.legend(handles=[Line2D([], [], ls="none", marker="s", ms=2.9, mfc="none", mec=OBS, mew=0.7, label="observed"),
                   Line2D([], [], color=BEG, marker="o", ms=2.7, mec="white", mew=0.3, lw=0.9, label="step 0"),
                   Line2D([], [], color=PRA, marker="o", ms=2.7, mec="white", mew=0.3, lw=0.9, label="step 1,000")],
          loc="upper left", fontsize=6.0, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.15, handlelength=1.3)

ax = axes[1]
f = task(ax, d["x_hyman"], "Hyman")
tb_axes(ax, xticks=[0, 1.5, 3], xlim=(-0.15, 3.15), yticks=[200, 450, 700], ylim=(150, 780))
ax.set_xlabel("information of the condition, bits", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)
slope_labels(ax, f, 2.05, 1.45, 30, 75, "{:.0f}")
ax.text(-0.05, 735, "24 condition means", fontsize=6.2, ha="left", va="center", color=GREY, **CLEAR)

ax = axes[2]
x = d["x_fitts"]
f = task(ax, x, "Fitts", use=x >= 2, obs=d["fitts_mt"], x_obs=d["fitts_id"])
tb_axes(ax, xticks=[1, 4, 7], xlim=(0.5, 7.5), yticks=[0, 400, 800], ylim=(0, 860))
ax.set_xlabel(r"$ID = \log_2(2D/W)$, bits", labelpad=1); ax.set_ylabel("movement time, ms", labelpad=1.5)
slope_labels(ax, f, 4.6, 3.9, 45, 110, "{:.1f}")
ax.annotate("no correction", xy=(1.08, 172), xytext=(1.5, 70), fontsize=6.2, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=RULE, shrinkA=1, shrinkB=1))

ax = axD
names = ("Hick", "Hyman", "Fitts")
for i, t in enumerate(names):
    r = d[f"ratio_{t}"]; lo, hi = np.percentile(r, [2.5, 97.5])
    ax.plot([i, i], [lo, hi], color=INK, lw=0.8, zorder=3)
    ax.plot(i, r.mean(), "o", ms=3.2, color=INK, mec="white", mew=0.3, zorder=4)
ax.axhline(H1 / H0, color=OBS, lw=0.9, zorder=2)
ax.axhline(L1 / L0, color=RULE, lw=0.7, ls=(0, (2.0, 1.6)), zorder=1)
ax.legend(handles=[Line2D([], [], color=OBS, lw=0.9, label=rf"$H$(1000) / $H$(0), {H1 / H0:.3f}"),
                   Line2D([], [], color=RULE, lw=0.7, ls=(0, (2.0, 1.6)), label=rf"$L$(1000) / $L$(0), {L1 / L0:.3f}")],
          loc="lower left", fontsize=6.0, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.15, handlelength=1.6)
ax.set_xlim(-0.5, 2.5); ax.set_xticks(range(3)); ax.set_xticklabels(names)
tb_axes(ax, yticks=[0.68, 0.695, 0.71], ylim=(0.672, 0.712), grid="y")
ax.set_ylabel("slope at step 1,000 / step 0", labelpad=1.5)

for ax_, L_, T_ in zip(axes + [axD], "ABCD", ("Hick", "Hyman", "Fitts", "One ratio")):
    S.head_over(ax_, L_, T_)
out = ROOT / "figures" / "figS2_static_cases.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, dict(zip("ABCD", ([a_] for a_ in axes + [axD]))), str(ROOT / "figures" / "panels" / "figS2"))
print({t: round(float(d[f"ratio_{t}"].mean()), 4) for t in names}, "H ratio", round(H1 / H0, 4), "L ratio", round(L1 / L0, 4))
print("written", out)
