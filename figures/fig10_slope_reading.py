"""Figure 10 of the manuscript: Leonard's practice read on the path of averaged information.

A, coupling of the model in May, taken as even over 20 elements, the case of Figure 2A. Above, the shares ranked from
the largest, read on the top axis, with the even share 1/N dashed. Below, the information blocks of width p_i and height
I(p_i) = -log2 p_i, whose area is H, with log2 N dashed. Between them one command, at the centre of the ranked shares,
fans out to the elements, each drawn above its block at the width of its share and reached by a branch whose width
follows the share (Assumption 1). B, the path of averaged information of the full model (fixed condition of Figure 5,
20 elements, mean and interquartile band of 300 runs), with Leonard's (1954) months placed by the ratio of each slope to
the slope of May (Equation 10): H = log2 20 x b(t)/b(May). The top axis names the months, and the right axis reads the
same path as Leonard's slope, beta H(t) with beta = b(May)/log2 20. The networks drawn into the plot, reached by dotted
arrows from May and December, are the coupling of the model at those two points, half transparent. The open circle
marks step 200, where selective elimination begins. The note gives the fall of the slope from May to December, which by
Equation 17 is also how much more time the six-choice task saved. The three markers at step 1,000 are the slopes of the
simulated performer of Supplemental Section S4.4 in the tasks of Hick, Hyman and Fitts, each divided by its slope at
step 0 and read on the left axis. C, as A for December, from the run whose path lies nearest the mean, at the step
where it reaches the averaged information of December, mirrored. An element whose block is taller than the row of
elements is drawn just above its block. Colour reads log2(p_i N) on the scale of Figure 2. Panels A and C and the
networks in B are the model at the averaged information that the data assign, not measured structures.
Inputs: data/transcribed/mowbray_rhoades1959_tables1_3.csv, code/sim_full.py and sim_out/sim_static_cases.npz (written
by code/sim_static_cases.py). The tables themselves are drawn in Supplemental Figure S3 (figures/figS3_practice_tables.py).
Run from the repository root: python figures/fig10_slope_reading.py sim_out   (about ten seconds)
"""
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(ROOT / "code"))
from style import WEIGHT, Sheet, INK, GREY, EDGE, RULE, LINK, TRACK, TONE, share_colour, save_panels, by_region  # noqa: E402
import sim_full as sf  # noqa: E402

logging.getLogger("matplotlib.axes._base").setLevel(logging.ERROR)   # round networks in panel B
SIM = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "sim_out"

# ---------------------------------------------------------------- data, Leonard (1954) as tabled by Mowbray and Rhoades
t13 = pd.read_csv(ROOT / "data" / "transcribed" / "mowbray_rhoades1959_tables1_3.csv")
MONTHS = ["May", "September", "December"]
lp = t13[t13.table == "I"].pivot(index="period_or_trial", columns="choices", values="mean_rt_ms").loc[MONTHS]
leo_b = (lp[6] - lp[3]).values                                      # ms per bit, log2 6 - log2 3 = 1
saved = (lp[3].iloc[0] - lp[3].iloc[-1], lp[6].iloc[0] - lp[6].iloc[-1])
assert list(leo_b) == [62, 51, 45] and saved == (12, 29)
r_leo = leo_b / leo_b[0]                                             # H(t)/H(May), Equation 10

# ---------------------------------------------------------------- the model path and the matched states
N, LAM, T_ELIM, T_END, RUNS = 20, 0.5, 200, 1000, 300
H0 = np.log2(N)
HT = np.array([sf.run(N=N, n_start=N, steps=T_END, lam=LAM, t_elim=T_ELIM, seed=s)[0] for s in range(RUNS)])
HM, Q1, Q3 = HT.mean(0), np.percentile(HT, 25, 0), np.percentile(HT, 75, 0)
target = H0 * r_leo                                                  # 4.322, 3.555, 3.137 bits
step = [0] + [int(np.argmax(HM <= h)) for h in target[1:]]
seed = int(np.argmin(((HT[:, :step[-1] + 1] - HM[:step[-1] + 1]) ** 2).mean(1)))   # the run nearest the mean path
step_run = [0] + [int(np.argmin(np.abs(HT[seed, max(0, s_ - 150): s_ + 151] - h))) + max(0, s_ - 150)
                  for s_, h in zip(step[1:], target[1:])]                 # where that run itself reaches each value
states = []
for s_ in step_run:
    if s_ == 0:
        A = np.zeros((N, N))
        for i in range(N):
            A[i, (i + 1) % N] += 1; A[(i + 1) % N, i] += 1
    else:
        A = sf.run(N=N, n_start=N, steps=s_, lam=LAM, t_elim=T_ELIM, seed=seed, return_adjacency=True)[3]
    states.append(A)

st = np.load(SIM / "sim_static_cases.npz")
RATIO = {t: float(st[f"ratio_{t}"].mean()) for t in ("Hick", "Hyman", "Fitts")}
BETA = leo_b[0] / H0                                                 # ms per command bit, Leonard in May


def shares(A):
    k = A.sum(1); p = k / k.sum(); I = np.full(N, np.nan); I[p > 0] = -np.log2(p[p > 0])
    return k, p, I, float(np.nansum(p * I))


# ---------------------------------------------------------------- layout
FW, FH = 6.5, 1.94
S = Sheet(FW, FH)
XL, XR, WS, XC, WC = 0.38, 4.70, 1.40, 2.30, 1.92                   # left, right and centre columns
YB, HB = 0.26, 0.74                                                 # spindles and information blocks
YA, HA = YB + HB, 0.47                                              # ranked shares, directly above
axA, axB = S.ax(XL, YA, WS, HA), S.ax(XL, YB, WS, HB)
axD, axE = S.ax(XR, YA, WS, HA), S.ax(XR, YB, WS, HB)
axC = S.ax(XC, YB, WC, HB + HA)                                     # B as high as A and C
FIRST, MID, LAST = TONE["small-world"], TONE["scale-free"], INK
PMAX = 0.30


def place(A, p):
    """Ring for even coupling; otherwise the fit hubs at the centre and each other element in the sector of its hub."""
    if np.allclose(p[p > 0], 1 / N):
        th = 2 * np.pi * np.arange(N) / N + np.pi / 2
        return {i: np.array([np.cos(th[i]), np.sin(th[i])]) for i in range(N)}
    hubs = list(np.argsort(p)[::-1][:3])
    ang = {h: np.pi / 2 + 2 * np.pi * j / 3 for j, h in enumerate(hubs)}
    pos = {h: 0.36 * np.array([np.cos(ang[h]), np.sin(ang[h])]) for h in hubs}
    group = {h: [] for h in hubs}
    for i in range(N):
        if i in hubs:
            continue
        group[hubs[int(np.argmax(A[i, hubs]))] if A[i, hubs].sum() > 0 else hubs[-1]].append(i)
    for h, mem in group.items():
        mem = sorted(mem, key=lambda i: -p[i])
        for j, i in enumerate(mem):
            a = ang[h] + (j - (len(mem) - 1) / 2) * (2 * np.pi / 3) / max(len(mem), 1) * 0.9
            pos[i] = np.array([np.cos(a), np.sin(a)])
    return pos


def draw_network(ax, A, p, size=1.0, alpha=1.0):
    pos = place(A, p)
    for i in range(N):
        for j in range(i + 1, N):
            if A[i, j] > 0:
                ax.plot(*zip(pos[i], pos[j]), color=LINK, lw=(0.25 + 0.30 * np.sqrt(A[i, j])) * size, alpha=0.85 * alpha,
                        solid_capstyle="round", zorder=2)
    for i in range(N):
        if p[i] == 0:                                               # no links, no share, as in Figure 2
            ax.scatter(*pos[i], s=size * 5, color="none", edgecolors=GREY, linewidths=0.45, alpha=alpha, zorder=5)
        else:
            ax.scatter(*pos[i], s=size * (4 + 520 * p[i]), color=[share_colour(np.log2(p[i] * N))], edgecolors=EDGE,
                       linewidths=0.4, alpha=alpha, zorder=5)


def ranked(ax, A, mirror):
    """Shares ranked from the top, as bars from the share axis read at the top, with the even share 1/N dashed."""
    k, p, I, H = shares(A)
    order = np.argsort(-p, kind="stable")
    for r, i in enumerate(order, 1):
        if p[i] > 0:
            ax.add_patch(Rectangle((0, r - 0.42), p[i], 0.84, fc=share_colour(np.log2(p[i] * N)), ec=EDGE, lw=0.2,
                                   zorder=3))
    ax.axvline(1 / N, color=INK, lw=0.7, ls=(0, (2.2, 1.6)), zorder=4)
    ax.set_xlim((PMAX, 0) if mirror else (0, PMAX)); ax.set_ylim(N + 0.6, 0.4)
    ax.xaxis.tick_top(); ax.xaxis.set_label_position("top")
    ax.set_xticks([0, 1 / N, 0.1, 0.2, 0.3]); ax.set_xticklabels(["0", "1/$N$", "0.1", "0.2", "0.3"])
    ax.tick_params(axis="x", length=1.8, pad=1, labelsize=6.2)
    ax.set_xlabel(r"share $p_i$", labelpad=1.5)
    ax.set_yticks([1, 20]); ax.set_ylabel("rank", labelpad=1.0)
    ax.yaxis.set_label_coords(1.19 if mirror else -0.20, 0.5)
    sd = "right" if mirror else "left"
    if mirror:
        ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
    for sp in ("left", "right"):
        ax.spines[sp].set_visible(sp == sd)
    ax.spines["top"].set_visible(True); ax.spines["top"].set_linewidth(0.6)
    ax.spines["bottom"].set_color(RULE); ax.spines["bottom"].set_linewidth(0.4)
    ax.patch.set_alpha(0.0)


def blocks(ax, A, mirror):
    """Information blocks of width p_i and height I(p_i), whose area is H, with each element above its block as a
    spindle of the same width."""
    k, p, I, H = shares(A)
    order = np.argsort(-p, kind="stable"); cum = 0.0; tips = []
    for i in order:
        if p[i] <= 0:
            continue
        c = share_colour(np.log2(p[i] * N)); xm = cum + p[i] / 2
        ax.add_patch(Rectangle((cum, 0), p[i], I[i], fc=c, ec=EDGE, lw=0.2, zorder=3))
        ys = max(YS, I[i] + 1.0)                                    # just above the H value, or above a tall block
        ax.add_patch(Ellipse((xm, ys), max(p[i] * 0.86, 0.006), 1.4, fc=c, ec=EDGE, lw=0.25, zorder=5))
        tips.append((xm, ys + 0.7, p[i])); cum += p[i]
    ax.plot([0, 1], [H0, H0], color=INK, lw=0.7, ls=(0, (2.2, 1.6)), zorder=4)
    ax.plot([0, 1], [H, H], color=INK, lw=1.1, zorder=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, YT); ax.set_yticks([0, H0, 9]); ax.set_yticklabels(["0", f"{H0:.2f}", "9"])
    ax.set_xticks([0, 0.5, 1])
    sd = "right" if mirror else "left"
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(sp == sd)
    ax.spines[sd].set_bounds(0, 9)
    ax.set_xticklabels(["0", "0.5", "1"]); ax.set_xlabel("cumulative share", labelpad=1)
    ax.set_ylabel(r"$I(p_i)$, bits", labelpad=1.5)
    if mirror:
        ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
    ax.yaxis.set_label_coords(1.19 if mirror else -0.20, 4.5 / YT)          # clear of the tick 4.32
    ax.text(0.98 if not mirror else 0.02, 4.75, f"$H$ = {H:.3f} bits", fontsize=6.0, ha="right" if not mirror else "left",
            va="bottom", zorder=6)
    ax.patch.set_alpha(0.0)
    return H, tips


def fan(x0, tips):
    """One command at the centre of the ranked shares, fanning out to the top of every spindle, line width by share."""
    ov = S.ax(x0, YB, WS, HB + HA); ov.set_xlim(0, WS); ov.set_ylim(0, HB + HA); ov.axis("off"); ov.set_zorder(0)
    cx, cy = WS / 2, HB + HA * 0.36
    for xm, yt, pi in tips:
        ov.plot([cx, xm * WS], [cy, yt / YT * HB], color=EDGE, lw=0.15 + 7.0 * pi, alpha=0.75, solid_capstyle="round")
    ov.plot(cx, cy, "o", ms=3.2, color=INK, zorder=6)
    ov.text(cx + 0.06, cy + 0.01, "command", fontsize=5.8, ha="left", va="center", zorder=6)


YS, YT = 7.95, 11.6
for a_, b_, x0, A_, mir in ((axA, axB, XL, states[0], False), (axD, axE, XR, states[2], True)):
    ranked(a_, A_, mir)
    H_, tips = blocks(b_, A_, mir)
    fan(x0, tips)
    a_.set_zorder(2); b_.set_zorder(2)

# ---------------------------------------------------------------- B, one path, two readings
ax = axC
X0, X1, Y0, Y1 = -40, 1100, 2.42, 4.62
for x_ in range(0, 1001, 200):
    ax.axvline(x_, color="#E9EDEF", lw=0.5, zorder=0.5)
for y_ in (2.5, 3.0, 3.5, 4.0, 4.5):
    ax.axhline(y_, color="#E9EDEF", lw=0.5, zorder=0.5)
tt = np.arange(T_END + 1)
ax.fill_between(tt, Q1, Q3, color=TRACK, lw=0, zorder=1)
ax.plot(tt, HM, color=EDGE, lw=1.2, zorder=2)
for m, h, s_, c in zip(MONTHS, target, step, (FIRST, MID, LAST)):
    ax.plot([s_, s_], [h, Y1], color=c, lw=0.7, ls=(0, (1.2, 1.2)), zorder=3)      # up to the month
    ax.plot([s_, X1], [h, h], color=c, lw=0.7, ls=(0, (1.2, 1.2)), zorder=3)      # across to the slope
    ax.plot(s_, h, "o", ms=5.2, color=c, mec="white", mew=0.7, zorder=6)
e = T_ELIM                                                           # where selective elimination begins
ax.plot(e, HM[e], "o", ms=4.0, color="white", mec=INK, mew=0.7, zorder=6)
ax.text(e - 25, HM[e] - 0.07, "elimination\nstarts", fontsize=5.8, ha="right", va="top", linespacing=0.95, zorder=7)
xa = 1072                                                            # May to December
ax.annotate("", xy=(xa, target[2]), xytext=(xa, target[0]), arrowprops=dict(arrowstyle="<->", lw=0.7, color=INK,
            mutation_scale=5, shrinkA=0, shrinkB=0), zorder=7)
ax.text(757, 3.94, rf"$\Delta H$ = {target[0] - target[2]:.2f} bits," + "\n" f"{leo_b[0] - leo_b[2]:.0f} ms per bit,"
        + "\n" f"{saved[1] - saved[0]:.0f} ms more" + "\n" "saved by six" + "\n" "choices", fontsize=5.4, ha="center", va="center",
        linespacing=1.05, zorder=7)
keys = []
for (t_, mk), x_ in zip((("Hick", "s"), ("Hyman", "D"), ("Fitts", "^")), (972, 1000, 1028)):
    ax.plot(x_, RATIO[t_] * H0, mk, ms=3.0, mfc="white", mec=INK, mew=0.6, zorder=7)
    keys.append(Line2D([], [], ls="none", marker=mk, ms=3.0, mfc="white", mec=INK, mew=0.6, label=t_))
ax.legend(handles=keys, title="simulated, step 1,000", loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=5.8,
          title_fontsize=5.8, handletextpad=0.15, borderaxespad=0.3, labelspacing=0.12, handlelength=0.9, ncol=3,
          columnspacing=0.45)
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
ax.set_xticks(range(0, 1001, 200)); ax.set_yticks([2.5, 3.0, 3.5, 4.0, 4.5])
for sp in ("left", "right"):
    ax.spines[sp].set_visible(False)
ax.spines["top"].set_visible(True)
ax.tick_params(axis="y", length=0, pad=2)
ax.set_xlabel("practice step $t$ of the model", labelpad=1); ax.set_ylabel("averaged information $H(t)$, bits", labelpad=1.5)
top = ax.secondary_xaxis("top"); top.set_xticks(step); top.set_xticklabels([m[:3] for m in MONTHS])
top.tick_params(length=1.8, pad=1, labelsize=6.4)
for lab, c in zip(top.get_xticklabels(), (FIRST, MID, LAST)):
    lab.set_color(EDGE if c == FIRST else c); lab.set_fontweight(WEIGHT)
top.spines["top"].set_linewidth(0.6)
sec = ax.secondary_yaxis("right", functions=(lambda h: h * BETA, lambda b: b / BETA))
sec.set_ticks(list(leo_b)); sec.tick_params(labelsize=6.4, length=0, pad=2)
for lab, c in zip(sec.get_yticklabels(), (FIRST, MID, LAST)):
    lab.set_color(EDGE if c == FIRST else c); lab.set_fontweight(WEIGHT)
sec.spines["right"].set_visible(False)
sec.set_ylabel(r"Leonard’s slope $\beta H(t)$, ms per bit", fontsize=6.8, labelpad=2)
NETS = (((205, 3.80, 420, 0.84), states[0], FIRST, (0, target[0]), (292, 4.24)),
        ((415, 2.40, 560, 0.90), states[2], LAST, (step[2], target[2]), (822, 3.00)))
for (x0, y0, w, h), A_, c, P, Q in NETS:                             # the coupling itself, drawn into the plot
    ins = ax.inset_axes([x0, y0, w, h], transform=ax.transData, zorder=1.5)
    ins.set_xlim(-1.2, 1.2); ins.set_ylim(-1.2, 1.2); ins.set_aspect("equal", adjustable="datalim"); ins.axis("off")
    ins.patch.set_alpha(0.0)
    draw_network(ins, A_, shares(A_)[1], size=0.78, alpha=0.55)
    ax.annotate("", xy=Q, xytext=P, arrowprops=dict(arrowstyle="-|>", lw=0.7, ls=(0, (1.2, 1.2)), color=c,
                mutation_scale=5, shrinkA=4, shrinkB=0), zorder=7)
S.head("A", "May, even coupling", XL + WS / 2, y=FH - 0.13)
S.head("B", "Leonard on the path of $H(t)$", XC + WC / 2, y=FH - 0.13)
S.head("C", "December, concentrated coupling", XR + WS / 2, y=FH - 0.13)
out = ROOT / "figures" / "fig10_slope_reading.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, by_region(S, {"A": (0, 0, XC - 0.12, FH), "B": (XC - 0.12, 0, XR - 0.06, FH), "C": (XR - 0.06, 0, FW, FH)}),
            str(ROOT / "figures" / "panels" / "fig10"))
print("Leonard slopes", leo_b.tolist(), "ratios", np.round(r_leo, 4).tolist(), "| H targets", np.round(target, 3).tolist())
print("time saved, three and six choices:", tuple(float(x) for x in saved), "| fall of the slope", leo_b[0] - leo_b[2], "ms per bit")
print("matched steps of the mean path", step, "| run shown, seed", seed, "at its own steps", step_run, "| its H there",
      [round(shares(A)[3], 3) for A in states])
print("model mean H at the matched steps", [round(float(HM[s_]), 3) for s_ in step], "| simulated performer ratios",
      {k: round(v, 4) for k, v in RATIO.items()}, "beta", round(BETA, 3))
print("written", out)
