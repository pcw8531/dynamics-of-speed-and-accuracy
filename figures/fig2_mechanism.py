"""Figure 2. The dynamic law as a scheme, T = a + beta H(t) log2 X, from even coupling (A) to hub coupling (D).

The classical laws (Equations 1 to 3) hold the coupling even, H = log2 N, and so describe A alone. The dynamic law
covers the whole path, because practice changes the coupling and with it H(t) and the slope b = beta H(t).

The figure has three layers. Above, the information of the elements: every element of a structure is a block of width
p_i (Equation 5, c_i = k_i) and height I(p_i) = -log2 p_i, so the total area of the blocks is H (Equations 6 and 7),
the solid line, with log2 N dashed. In the middle, the structures themselves, with node area by share and link width by
the coupling of the pair, each joined by dotted arrows to two of its blocks. Below, the law with its terms named,
between the selection of a command at A and at D. The words classical laws, static stand inside the chain of A, and a
dotted bar over the law, open towards it, names the dynamic law.

One diverging scale reads the share of every element, on log2(p_i N): blue below an even share, pale yellow at the
even share 1/N, and yellow to red above it, so hubs are the reddest elements. It is the palette of Park (2026), the
one of figures/style.py, so every figure that follows Figure 2 has the same colours.

A, the novice, even coupling: the chain of 20 elements and 20 links (H = log2 20). D, the expert, hub coupling,
mirrored: the scale-free structure of 20 elements. Beside each of the two structures stands its H in bits, and both
are drawn large with nodes that let through what lies under them. Between A and D, unlabelled, the small-world and
random structures.
B, H(t) under the practice rules of the full model (code/sim_full.py), whose starting ring is the chain of A: formation
by eta_i k_i, selective elimination by 1/eta per link end from step 200, lambda = 0.5 with g = 1 on three elements; 300 runs,
mean and interquartile band. C, the law at A, at D and at the end of practice (Equation 10).
Below A and below D, the selection of a command as a race (section From Averaged Information to the Slope;
Supplemental Material, Section S8), on one scale: mean paths from the head start h_i = ln p_i to the criterion theta,
at the rate v, so an element is selected after (theta - h_i) / v. Left, every element of A. Right, mirrored, the hub
and a leaf of D. Line width reads the share, the probability that a command addresses the element, and the tick marked
A is the time at which every element of A is selected. Drawn with theta = 3; the times are printed. The head-start
axis is in natural units, with ticks at the criterion 3, at 0 and at -3, the head start ln(1/20) of every element of A.

The four structures of 20 elements and 20 links are the sparse pair of the deposit of Park (2026) with the two
structures between them, all made by generators.make of that deposit at mean degree 2:
  regular       ring lattice, a closed chain (H = log2 N)
  small-world   the same generator at rewiring probability 0.1 (Watts & Strogatz, 1998)
  random        20 links placed at random (Erdos & Renyi, 1959)
  scale-free    preferential attachment with m = 1 from a ring seed of 4 elements (Barabasi & Albert, 1999)
For each generated structure the realisation shown is the one nearest the mean H of seeds 42 to 141. The chain is
drawn along a smooth closed contour, and the small-world and random structures on the same places, so the loss of
the chain's order shows; the scale-free one is drawn in rings around its strongest element, as in the earlier figure.
Every drawing is stretched sideways by 1.2 with round nodes. An element left without links is drawn as an open point.
Writes figures/fig2_mechanism.png, each lettered panel and the centre strip to figures/panels, and prints every
number the figure shows.
Usage: python -I figures/fig2_mechanism.py <path to the deposit clone at commit 0660a649f5e0>
"""
import logging, os, subprocess, sys
from collections import Counter
from pathlib import Path
import numpy as np
import networkx as nx
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.colorbar  # noqa: E402,F401
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle, FancyArrowPatch, ConnectionPatch
from matplotlib.lines import Line2D
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
ROOT = Path(__file__).resolve().parent.parent
DEP = Path(sys.argv[1]).resolve()
head_ = subprocess.run(["git", "-C", str(DEP), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
assert head_.startswith("0660a649f5e0"), f"deposit is at {head_}"
sys.path.insert(0, str(DEP / "code")); sys.path.insert(0, str(ROOT / "code"))
from generators import make, SEEDS                                   # noqa: E402  Park (2026)
import sim_full as sf                                                # noqa: E402  the full model

mpl.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 7.5,
                     "mathtext.fontset": "stix", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 1.8, "ytick.major.size": 1.8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 7.5, "savefig.dpi": 600})
# regular type and narrow trend lines, as in figures/style.py (FIG_WEIGHT=bold gives bold type)
BOLD = os.environ.get("FIG_WEIGHT", "plain") == "bold"
WEIGHT, _W = ("bold", ":bold") if BOLD else ("normal", "")
mpl.rcParams.update({"font.weight": WEIGHT, "axes.labelweight": WEIGHT, "axes.titleweight": WEIGHT,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Serif" + _W,
                     "mathtext.it": "Liberation Serif:italic" + _W, "mathtext.bf": "Liberation Serif:bold",
                     "mathtext.sf": "Liberation Serif" + _W, "mathtext.fallback": "stix"})
# One diverging scale for the share of every element, read on log2(p_i N) in bits on +-2.716: blue below an even share,
# pale yellow at exactly 1/N, yellow to red above it. This is the palette of Park (2026), as in figures/style.py.
INK, GREY, GUIDE, RULE, TRACK = "#1A1A1A", "#8A8A8A", "#455A64", "#8C989E", "#E3E7E9"
LINK, EDGE = "#C3CBD0", "#37474F"                                     # links; outlines, stems and arrows
PALE = TRACK
ENS_TONE = {"regular": "#B9C4C9", "small-world": "#93A8B1", "random": "#6E8C99", "scale-free": "#48707F"}
VLIM = 2.7157337518864986
_AT = lambda b: (b + VLIM) / (2 * VLIM)                               # place of log2(p_i N) = b on the scale
_HOT = plt.get_cmap("hot_r")(np.linspace(0.06, 0.94, 256))
_COOL = plt.get_cmap("Blues")(np.linspace(0.92, 0.06, 256))
CMAP = LinearSegmentedColormap.from_list("hot_div", np.vstack([_COOL, _HOT]))   # the palette of Park (2026), as in figures/style.py

FW, FH = 6.5, 3.05
Y_TITLE = FH - 0.14


# ---------------------------------------------------------------- the four structures
N, KMEAN = 20, 2
KINDS = ("regular", "small-world", "random", "scale-free")


def readings(G):
    k = np.array([G.degree(n) for n in sorted(G.nodes())], float); p = k / k.sum()   # Equation 5, c_i = k_i
    I = -np.log2(np.where(p > 0, p, np.nan))                                         # information of element i
    return dict(k=k, p=p, I=I, H=float(np.nansum(p * I)))                             # Equation 6


NET, ENS = {"regular": make("regular", N, KMEAN, SEEDS[0])}, {}
for kind in KINDS[1:]:
    Gs = [make(kind, N, KMEAN, s) for s in SEEDS]
    ENS[kind] = float(np.mean([readings(G)["H"] for G in Gs]))
    NET[kind] = min(Gs, key=lambda G: abs(readings(G)["H"] - ENS[kind]))
assert all(G.number_of_edges() == N * KMEAN // 2 for G in NET.values())
R = {k: readings(G) for k, G in NET.items()}
CEIL = np.log2(N)
A_MS, BETA = 12.8, 94.7 / CEIL                                       # the chain reproduces Fitts' slope
assert abs(R["regular"]["H"] - CEIL) < 1e-12

# ---------------------------------------------------------------- practice from A: the fixed regime of the full model
LAM, T_ELIM, T_END, RUNS = 0.5, 200, 1000, 300
HT_ = np.array([sf.run(N=N, n_start=N, steps=T_END, lam=LAM, t_elim=T_ELIM, seed=s)[0] for s in range(RUNS)])
HMEAN = HT_.mean(0)
PASS = {k: int(np.argmax(HMEAN <= R[k]["H"])) for k in KINDS[1:]}

PHI = max(R[k]["p"].max() for k in KINDS)


def pcol(I, n):
    """Colour of an element of information I among n, by log2(p_i n) = log2 n - I(p_i), the bits it lies above the even share."""
    b = float(np.clip(np.log2(n) - I, -VLIM, VLIM))
    return CMAP(_AT(0.0 if abs(b) < 1e-9 else b))


# ---------------------------------------------------------------- layouts
def contour(n, seed=5):
    """A smooth closed contour through the chain, so it reads as a chain and not as a circle."""
    rng = np.random.default_rng(seed); th = 2 * np.pi * np.arange(n) / n
    r = 1 + sum(a * np.cos(m * th + ph) for m, a, ph in zip((2, 3, 4), (0.26, 0.18, 0.08), rng.uniform(0, 2 * np.pi, 3)))
    xy = np.c_[r * np.cos(th), r * np.sin(th)]; xy -= xy.mean(0); xy /= np.abs(xy).max()
    return {i: xy[i] * 0.95 for i in range(n)}


def radial(G, w, gap=225.0):
    """Rings by distance from the element with the largest share, each subtree given an angular span in proportion
    to its size (the layout of the earlier figure), turned so that a gap between subtrees faces `gap` degrees."""
    root = max(G.nodes(), key=lambda z: (w[z], G.degree(z), -z))
    parent, order, depth, seen, queue = {root: None}, [root], {root: 0}, {root}, [root]
    while queue:
        u = queue.pop(0)
        for v in sorted(G.neighbors(u), key=lambda z: (-w[z], -G.degree(z), z)):
            if v not in seen:
                seen.add(v); parent[v] = u; depth[v] = depth[u] + 1; order.append(v); queue.append(v)
    kids = {n: [] for n in order}
    for n in order[1:]:
        kids[parent[n]].append(n)
    size = {}
    for n in reversed(order):
        size[n] = max(1, sum(size[c] for c in kids[n])) if kids[n] else 1
    span = {root: (0.0, 2 * np.pi)}
    for n in order:
        a0, a1 = span[n]; tot = sum(size[c] for c in kids[n]) or 1; acc = a0
        for c in kids[n]:
            wd = (a1 - a0) * size[c] / tot; span[c] = (acc, acc + wd); acc += wd
    dmax = max(depth.values()) or 1
    ang = {n: 0.5 * sum(span[n]) for n in order}; rad = {n: np.sqrt(depth[n] / dmax) for n in order}
    first = sorted(ang[c] for c in kids[root])                       # the widest gap between first-ring subtrees
    gaps = [((first[(i + 1) % len(first)] - first[i]) % (2 * np.pi), first[i]) for i in range(len(first))]
    wid, start = max(gaps); turn = np.deg2rad(gap) - (start + wid / 2)
    return {n: np.array([rad[n] * np.cos(ang[n] + turn), rad[n] * np.sin(ang[n] + turn)]) * 0.95 for n in order}


LOOP = contour(N)
POS = {"regular": LOOP, "small-world": LOOP, "random": LOOP,
       "scale-free": radial(NET["scale-free"], R["scale-free"]["p"])}


def ax_in(x, y, w, h):
    return fig.add_axes([x / FW, y / FH, w / FW, h / FH])


def head(letter, title, x_centre, y=Y_TITLE, gap=0.06):
    """Panel letter set right before its title, the two centred together over the panel."""
    rend = fig.canvas.get_renderer()
    widths = []
    for txt, kw in ((letter, dict(fontsize=9, fontweight="bold")), (title, dict(fontsize=8))):
        t_ = fig.text(0, 0, txt, **kw); widths.append(t_.get_window_extent(rend).width / fig.dpi); t_.remove()
    x0 = x_centre - (widths[0] + gap + widths[1]) / 2
    fig.text(x0 / FW, y / FH, letter, fontsize=9, fontweight="bold", ha="left", va="baseline")
    fig.text((x0 + widths[0] + gap) / FW, y / FH, title, fontsize=8, ha="left", va="baseline")


def heading(letter, title, x, outer, gap=0.06):
    """Panel letter set right before its title, the pair starting at the left edge (outer = "left") or ending at the
    right edge (outer = "right") of the figure."""
    rend = fig.canvas.get_renderer(); y = FH - 0.14
    widths = []
    for txt, kw in ((letter, dict(fontsize=9, fontweight="bold")), (title, dict(fontsize=8))):
        t_ = fig.text(0, 0, txt, **kw); widths.append(t_.get_window_extent(rend).width / fig.dpi); t_.remove()
    x0 = x if outer == "left" else x - (widths[0] + gap + widths[1])
    fig.text(x0 / FW, y / FH, letter, fontsize=9, fontweight="bold", ha="left", va="baseline")
    fig.text((x0 + widths[0] + gap) / FW, y / FH, title, fontsize=8, ha="left", va="baseline")


def link_width(pu, pv, size):
    return (0.25 + 1.3 * (np.sqrt(pu * pv) / PHI) ** 0.45) * np.sqrt(size)


def wide(ax, pos, sx, lim):
    """Places stretched sideways by sx in an axes wider than tall, drawn at equal scale so nodes stay round."""
    bb = ax.get_position(); ar = (bb.width * FW) / (bb.height * FH)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-lim * ar, lim * ar); ax.set_ylim(-lim, lim)
    return {n: np.array([xy[0] * sx, xy[1]]) for n, xy in pos.items()}


def network(ax, kind, size=1.0, sx=1.0, alpha=1.0):
    """Node area by share p_i, link width by the coupling of the pair, node shade by information I(p_i).
    With alpha below 1 the fill of a node lets through what lies under it, and its outline stays solid."""
    G, Rk = NET[kind], R[kind]
    pos = wide(ax, POS[kind], sx, 1.1)
    for u, v in G.edges():
        ax.plot(*zip(pos[u], pos[v]), color=LINK, lw=link_width(Rk["p"][u], Rk["p"][v], size), solid_capstyle="round", zorder=2)
    for n in range(N):
        if Rk["p"][n] == 0:                                          # no links, no share
            ax.scatter(*pos[n], s=size * 6, color="white", edgecolors=GREY, linewidths=0.45, zorder=5)
            continue
        ax.scatter(*pos[n], s=size * (3 + 560 * Rk["p"][n]), color=[pcol(Rk["I"][n], N)[:3] + (alpha,)], edgecolors=EDGE,
                   linewidths=0.45, zorder=5)
    return pos


YMAX = 6.0


def guide(axA, xyA, axB, xyB, rad=0.0, shrinkB=1, style=None):
    fig.add_artist(ConnectionPatch(xyA=xyA, coordsA=axA.transData, xyB=xyB, coordsB=axB.transData, arrowstyle="-|>",
                                   mutation_scale=6, lw=0.7, color=GUIDE, ls=(0, (2.0, 1.6)), zorder=9,
                                   shrinkA=4, shrinkB=shrinkB, connectionstyle=style or f"arc3,rad={rad}"))


def fig_arrow(x0, x1, y):
    fig.add_artist(FancyArrowPatch((x0 / FW, y / FH), (x1 / FW, y / FH), transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=8, lw=0.9, color=EDGE, ls=(0, (2.4, 1.8)), zorder=9))


def top_bottom_box(ax):
    for sp in ("left", "right"):
        ax.spines[sp].set_visible(False)


def scope(x0, x1, y, label, tick=0.07, at=None, fontsize=7.2, ls="-"):
    """A bar beside what it covers, its ends turned towards it by `tick`, its label set in the bar, at its middle or at `at`."""
    fig.add_artist(Line2D(np.array([x0, x0, x1, x1]) / FW, np.array([y + tick, y, y, y + tick]) / FH,
                          transform=fig.transFigure, color=INK, lw=0.7, ls=ls))
    fig.text(((x0 + x1) / 2 if at is None else at) / FW, y / FH, label, fontsize=fontsize, ha="center", va="center",
             bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none"))


def ftext(x, y, s, **kw):
    return fig.text(x / FW, y / FH, s, **kw)




def info_blocks(ax, Rk, n, mirror=False, order=None, ylabel=True):
    """Blocks of width p_i and height I(p_i), in order of rising share from the outer edge inward (ties in the
    given order). The right-hand panel is mirrored, so its shares run from right to left. Returns the foot of each
    block, where the arrow from its element arrives."""
    cum, foot = 0.0, {}
    for i in (np.argsort(Rk["p"], kind="stable") if order is None else order):
        ax.add_patch(Rectangle((cum, 0), Rk["p"][i], Rk["I"][i], facecolor=pcol(Rk["I"][i], n), edgecolor=EDGE, lw=0.25, zorder=3))
        foot[i] = (cum + Rk["p"][i] / 2, 0.0); cum += Rk["p"][i]
    ax.plot([0, 1], [np.log2(n)] * 2, color=RULE, lw=0.7, ls=(0, (2.4, 1.6)), zorder=6)
    ax.plot([0, 1], [Rk["H"]] * 2, color=INK, lw=1.3, zorder=8)
    ax.set_xlim((1, 0) if mirror else (0, 1)); ax.set_ylim(0, YMAX); ax.set_yticks([0, 2, 4])
    ax.set_xticks([0, 1]); ax.set_xticklabels(["0", "1"])
    ax.tick_params(axis="x", pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["left" if mirror else "right"].set_visible(False)
    if mirror:
        ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
    if ylabel:
        ax.set_ylabel("information\n" + r"$I_i$, bits", labelpad=1.5, fontsize=7.0, linespacing=0.95)
    ax.text(0.5, -0.075, "share of coupling", transform=ax.transAxes, fontsize=7.5, ha="center", va="top")
    return foot


# ---------------------------------------------------------------- rows and columns
Y1, H1 = 2.04, 0.72                                                   # above: the blocks of A and D, and B and C
Y2, H2 = 1.06, 0.74                                                   # middle: the structures
W_P = 1.60
X_A, X_D = 0.44, FW - 0.44 - W_P
XC_A, XC_D = X_A + W_P / 2, X_D + W_P / 2
XM = FW / 2
SX = 1.2                                                              # sideways stretch of every network drawing
S_N, A_N = 1.45, 0.85                                                 # the structures of A and D: scale, and opacity of the nodes
W_N, H_N, Y_N = 1.10 * S_N, H2 * S_N, Y2 + H2 / 2 - H2 * S_N / 2      # drawn large, over the edges of the rows beside them

fig = plt.figure(figsize=(FW, FH))

# ---------------------------------------------------------------- A, the novice, even coupling: the blocks above, the chain below
heading("A", "Novice, even coupling", 0.04, "left")
RA = R["regular"]
by_x = sorted(range(N), key=lambda n: LOOP[n][0])                    # equal shares, ties ordered by position
pA = ax_in(X_A, Y1, W_P, H1); footA = info_blocks(pA, RA, N, order=by_x)
pA.text(0.5, CEIL + 0.3, r"$H = \log_2 N$", fontsize=7.0, ha="center", va="bottom")
aA = ax_in(XC_A - W_N / 2, Y_N, W_N, H_N); aA.set_zorder(4); pa = network(aA, "regular", size=S_N ** 2, sx=SX, alpha=A_N)
up = sorted((n for n in range(N) if LOOP[n][1] > 0.25), key=lambda n: LOOP[n][0])   # two upper elements, left and right
gA = (up[0], up[-1])
for n in gA:
    guide(aA, pa[n], pA, footA[n])

# ---------------------------------------------------------------- D, the expert, hub coupling, mirrored: the blocks above, the structure below
heading("D", "Expert, hub coupling", FW - 0.04, "right")
RD = R["scale-free"]; pos_d = POS["scale-free"]
hub = int(np.argmax(RD["p"]))
leaves = np.flatnonzero(RD["k"] == RD["k"].min())
weak = int(max(leaves, key=lambda n: pos_d[n][0] + pos_d[n][1]))     # a leaf toward the upper right
pD = ax_in(X_D, Y1, W_P, H1)
footD = info_blocks(pD, RD, N, mirror=True, order=sorted(range(N), key=lambda n: (RD["p"][n], n != weak)))
pD.text(0.62, CEIL + 0.3, r"$H < \log_2 N$", fontsize=7.0, ha="center", va="bottom")
aD = ax_in(XC_D - W_N / 2, Y_N, W_N, H_N); aD.set_zorder(4); pd_ = network(aD, "scale-free", size=S_N ** 2, sx=SX, alpha=A_N)
guide(aD, pd_[hub], pD, footD[hub], rad=-0.12)
guide(aD, pd_[weak], pD, footD[weak])

# the one value of each structure, beside it on its outer side: the bits per command
X_VAL = 0.34
for x_, R_ in ((X_VAL, RA), (FW - X_VAL, RD)):
    ftext(x_, Y2 + H2 / 2, rf"$H$ = {R_['H']:.2f}" + "\nbits", fontsize=7.6, ha="center", va="center", linespacing=1.05)

# ---------------------------------------------------------------- B, H(t) from A, with the levels of A and D
X_B, W_B, X_C, W_C = 2.50, 0.52, 3.22, 0.84
head("B", "Practice", X_B + W_B / 2)
bB = ax_in(X_B, Y1, W_B, H1)
t = np.arange(T_END + 1)
bB.fill_between(t, np.percentile(HT_, 25, axis=0), np.percentile(HT_, 75, axis=0), color=PALE, lw=0, zorder=1)
bB.plot(t, HMEAN, color=INK, lw=1.4, zorder=4)
for k, lab in (("regular", "A"), ("scale-free", "D")):
    bB.axhline(R[k]["H"], color=GUIDE, lw=0.8, ls=(0, (1.2, 1.4)), zorder=2)
    bB.text(T_END * 1.04, R[k]["H"], lab, fontsize=7, fontweight=WEIGHT, ha="left", va="center", clip_on=False)
for k in ("small-world", "random", "scale-free"):                    # where practice passes each structure
    bB.plot([PASS[k]], [R[k]["H"]], "o", ms=2.6, color=ENS_TONE[k], mec=EDGE, mew=0.5, zorder=6)
bB.set_xlim(0, T_END); bB.set_ylim(2.85, 4.48); bB.set_yticks([3.0, 4.0]); bB.set_yticklabels(["3", "4"])
bB.set_xticks([0, 1000]); bB.set_xticklabels(["0", "1,000"])
top_bottom_box(bB)
bB.set_xlabel("practice step $t$", labelpad=0.8); bB.set_ylabel(r"$H(t)$, bits", labelpad=1.5)

# ---------------------------------------------------------------- C, the law at A, at D and at the end of practice
head("C", r"Slope, $b = \beta H(t)$", X_C + W_C / 2)
cC = ax_in(X_C, Y1, W_C, H1)
xb = np.array([0.0, 4.0])
rows = [(CEIL, dict(color=ENS_TONE["regular"], lw=1.3), "A"), (RD["H"], dict(color=ENS_TONE["scale-free"], lw=1.3), "D"),
        (HMEAN[-1], dict(color=INK, lw=1.6), "$t$ = 1,000")]
for H_, kw, lab in rows:
    cC.plot(xb, A_MS + BETA * H_ * xb, label=lab, **kw)
cC.annotate("", xy=(4.0, A_MS + BETA * HMEAN[-1] * 4.0), xytext=(4.0, A_MS + BETA * CEIL * 4.0),
            arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.7, color=INK), annotation_clip=False)
cC.legend(loc="upper left", fontsize=6.0, frameon=False, handlelength=0.9, handletextpad=0.3, borderaxespad=0.0,
          labelspacing=0.0)
cC.set_xlim(0, 4.0); cC.set_xticks([0, 4]); cC.set_ylim(0, 480); cC.set_yticks([0, 200, 400])
top_bottom_box(cC); cC.yaxis.tick_right(); cC.yaxis.set_label_position("right")
cC.set_xlabel(r"$\log_2 X$, bits", labelpad=0.8); cC.set_ylabel(r"$T$, ms", labelpad=1.5)

# ---------------------------------------------------------------- between A and D, small, the small-world and random structures
S_M = 0.42
Y_M = Y2 + H2 / 2 - S_M / 2
sz = (S_M / H2) ** 2
aS = ax_in(XM - 0.55 - S_M / 2, Y_M, S_M, S_M); network(aS, "small-world", size=sz)
aR = ax_in(XM + 0.55 - S_M / 2, Y_M, S_M, S_M); network(aR, "random", size=sz)
y_mid = Y_M + S_M / 2
ftext(XM, y_mid, r"$\cdots$", fontsize=12, ha="center", va="center")
fig_arrow(XC_A + 0.76, XM - 0.55 - S_M / 2 - 0.02, y_mid); fig_arrow(XM + 0.55 + S_M / 2 + 0.02, XC_D - 0.76, y_mid)
ftext(XM, y_mid - 0.10, "practice", fontsize=6.8, ha="center", va="top")

# ---------------------------------------------------------------- the law in the middle, its terms named
y_eq = 0.71
fig.canvas.draw(); rend = fig.canvas.get_renderer()
parts = [(r"$T\ =\ a\ +\ $", None, 0), (r"$\beta$", "ms per bit", -0.40), (r"$\,H(t)$", "bits per command", 0.0),
         (r"$\ \log_2 X$", "steps of the act", 0.42)]
xs, x0 = [], 0.0
for s_, _, _ in parts:
    t_ = fig.text(0, 0, s_, fontsize=12); w_ = t_.get_window_extent(rend).width / fig.dpi; t_.remove(); xs.append((x0, w_)); x0 += w_
left = XM - x0 / 2
for (s_, lab, dx), (xo, w_) in zip(parts, xs):
    ftext(left + xo, y_eq, s_, fontsize=12, ha="left", va="center")
    if lab:
        xc = left + xo + w_ / 2
        fig.add_artist(FancyArrowPatch(((xc + dx) / FW, (y_eq - 0.24) / FH), (xc / FW, (y_eq - 0.12) / FH),
                                       transform=fig.transFigure, arrowstyle="-|>", mutation_scale=5, lw=0.5, color=GREY))
        ftext(xc + dx, y_eq - 0.26, lab, fontsize=6.4, ha="center", va="top")

# ---------------------------------------------------------------- below A and below D, the selection of a command as a race
# Each element races from its head start h_i = ln p_i to the criterion theta at the rate v, so it is selected after
# (theta - h_i) / v. Mean paths on one scale in both panels, the width of each path by the share of its element.
THETA, UNIT = 3.0, BETA / np.log(2.0)                                # criterion, and ms per unit of 1 / v
T_MAX = 252.0
Y_E, H_E, W_E = 0.38, 0.47, 1.14
LO = np.log(RD["p"][weak]) - 0.75
T_ARR = {"A": (THETA - np.log(1.0 / N)) * UNIT}


def race_panel(x0, marks, mirror=False):
    """Mean paths of the elements in marks, a list of (share, information, label, place of the label along the path,
    positive for a label above the path, negative for one below it, and None for a label at its tick on the time axis).
    The right-hand panel is mirrored, as the blocks of D are, so its time runs from right to left."""
    ax = ax_in(x0, Y_E, W_E, H_E)
    ax.axhline(THETA, color=INK, lw=0.7, ls=(0, (2.4, 1.6)), zorder=2)
    ax.text(T_MAX, THETA + 0.3, r"criterion $\theta$", fontsize=6.2, ha="left" if mirror else "right", va="bottom")
    for p_, I_, lab, f_ in marks:
        h_ = np.log(p_); T_ARR[lab] = (THETA - h_) * UNIT
        ax.plot([0, T_ARR[lab]], [h_, THETA], color=GUIDE, lw=0.45 + 9.0 * p_, zorder=3, solid_capstyle="round")
        ax.scatter([0], [h_], s=0.55 * (3 + 560 * p_), color=[pcol(I_, N)], edgecolors=EDGE, linewidths=0.45, zorder=5, clip_on=False)
        ax.plot([T_ARR[lab]] * 2, [LO, THETA], color=RULE, lw=0.5, ls=(0, (1.0, 1.4)), zorder=1)
        if f_ is not None:
            ax.text(abs(f_) * T_ARR[lab], h_ + abs(f_) * (THETA - h_) + (0.8 if f_ > 0 else -0.55), lab, fontsize=6.2,
                    color=EDGE, ha="right" if f_ > 0 else "left", va="bottom" if f_ > 0 else "top")
    ax.set_xlim((T_MAX, 0) if mirror else (0, T_MAX)); ax.set_ylim(LO, THETA + 0.2)
    at_tick = sorted([("A", T_ARR["A"])] + [(lab, T_ARR[lab]) for _, _, lab, f_ in marks if f_ is None], key=lambda z: z[1])
    ax.set_xticks([v for _, v in at_tick]); ax.set_yticks([-3, 0, 3])   # natural units: the criterion at 3, an even share of 1/20 at -3.0
    if mirror:
        ax.yaxis.tick_right()
    for tl, (lab, _) in zip(ax.set_xticklabels([lab for lab, _ in at_tick]), at_tick):
        if lab == "leaf":                                            # the last label set beyond its tick
            tl.set_ha("right" if mirror else "left")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=6.4, pad=1.0, length=2.4)
    ax.tick_params(axis="y", labelsize=6.4, pad=3.2 if mirror else 2.5, length=2.0)   # clear of the start markers
    ax.set_xlabel("time to select", labelpad=0.5, fontsize=6.8)
    ax.text(1.375 if mirror else -0.21, 0.42, r"head start $\ln p_i$", transform=ax.transAxes, fontsize=6.2, rotation=90,
            ha="center", va="center")                                   # on the outer side of each panel, as for the blocks
    return ax


X_ED = X_D + 0.34
eA = race_panel(FW - X_ED - W_E, [(1.0 / N, CEIL, "every element", -0.36)])
eD = race_panel(X_ED, [(RD["p"][hub], RD["I"][hub], "hub", None), (RD["p"][weak], RD["I"][weak], "leaf", None)], mirror=True)
# each structure joined to its race, as it is to its blocks: an element of A, and the hub and a leaf of D
low_a = min(range(N), key=lambda n: pa[n][0] + pa[n][1])
guide(aA, pa[low_a], eA, (0, np.log(1.0 / N)), rad=-0.10, shrinkB=3.5)
guide(aD, pd_[hub], eD, (0, np.log(RD["p"][hub])), rad=0.10, shrinkB=5.5)
# the leaf starts under the hub, so its arrow is led round to arrive from the outer side
low_leaf = int(max(leaves, key=lambda n: pd_[n][0] - 0.5 * pd_[n][1]))   # a leaf at the lower right
fig.canvas.draw()
inch = lambda ax_, xy: fig.dpi_scale_trans.inverted().transform(ax_.transData.transform(xy))
a_, b_ = inch(aD, pd_[low_leaf]), inch(eD, (0, np.log(RD["p"][low_leaf])))
from matplotlib.path import Path as MPath  # noqa: E402
bend = MPath(np.array([a_ + np.array([0.012, -0.045]), np.array([b_[0] + 0.33, b_[1] + 0.58]),
                       np.array([b_[0] + 0.36, b_[1] - 0.22]), b_ + np.array([0.026, -0.036])]) / np.array([FW, FH]),
             [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])   # in parts of the figure, so a cut panel keeps it
fig.add_artist(FancyArrowPatch(path=bend, transform=fig.transFigure, arrowstyle="-|>", mutation_scale=6, lw=0.7,
                               color=GUIDE, ls=(0, (2.0, 1.6)), zorder=9, shrinkA=0, shrinkB=0))

# ---------------------------------------------------------------- what each law covers: the classical laws A alone, the dynamic law A to D
# the classical laws cover A alone: named inside the chain, at the centre of the area it encloses
ring = np.array([inch(aA, pa[n]) for n in range(N)]); x_, y_ = ring[:, 0], ring[:, 1]
cr = x_ * np.roll(y_, -1) - np.roll(x_, -1) * y_
cx, cy = ((x_ + np.roll(x_, -1)) * cr).sum() / (3 * cr.sum()), ((y_ + np.roll(y_, -1)) * cr).sum() / (3 * cr.sum())
ftext(cx, cy, "classical laws,\nstatic", fontsize=7.2, ha="center", va="center", linespacing=1.05, zorder=10)
# the dynamic law covers the whole path: named over the law, on a dotted bar that opens towards it
scope(2.02, FW - 2.02, 1.00, r"dynamic law, $H(t)$ from A to D", tick=-0.07, fontsize=9.0, ls=(0, (2.0, 1.6)))

# the trend lines of B and C set narrow, as in figures/style.py (bolden): lines of at least 0.95 pt in axes with ticks
# ---------------------------------------------------------------- the finished figure, and each panel on its own
FACTOR = 1.4 if BOLD else 0.75
for ax_ in (bB, cC):
    for ln in ax_.get_lines():
        if ln.get_linestyle() not in ("None", "", " ") and ln.get_linewidth() >= 0.95:
            ln.set_linewidth(ln.get_linewidth() * FACTOR)
    leg = ax_.get_legend()
    for h_ in (leg.legend_handles if leg is not None else []):
        if hasattr(h_, "set_linewidth") and h_.get_linestyle() not in ("None", "", " ") and h_.get_linewidth() >= 0.95:
            h_.set_linewidth(h_.get_linewidth() * FACTOR)
OUT = Path(os.environ.get("FIG2_OUT", ROOT / "figures" / "fig2_mechanism.png"))
fig.savefig(OUT, dpi=600)

# each lettered panel on its own: the artists of the figure are sorted by where their centres lie, the others are hidden
# while the panel is saved, and the crop is the union of what remains. A and D take their column whole, the blocks, the
# structure and the selection of a command. The centre strip, practice between A and D and the law under its dotted bar,
# has no letter and is written as fig2_centre.
from matplotlib.transforms import Bbox  # noqa: E402
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.dpi_scale_trans.inverted()
X_END_A, X_BEG_D, Y_LOW = X_A + W_P + 0.10, X_D - 0.10, Y1 - 0.42


def panel_of(bb):
    xc, yc = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
    if yc < 0.18:
        return None
    if xc < X_END_A:
        return "A"
    if xc > X_BEG_D:
        return "D"
    if yc > Y_LOW:
        return "B" if xc < X_C - 0.02 else "C"
    return "centre"


def extent(a_):
    bb = a_.get_tightbbox(rend) if hasattr(a_, "get_tightbbox") else None
    return (bb if bb is not None else a_.get_window_extent(rend)).transformed(inv)


arts = [a_ for a_ in list(fig.axes) + list(fig.texts) + list(fig.artists) + list(fig.lines) + list(fig.patches)
        if a_.get_visible()]
boxes = [extent(a_) for a_ in arts]
where = [panel_of(bb) for bb in boxes]
for name in ("A", "B", "C", "D", "centre"):
    keep = [i for i, w in enumerate(where) if w == name]
    for i, a_ in enumerate(arts):
        a_.set_visible(i in keep)
    bb = Bbox.union([boxes[i] for i in keep]).padded(0.04)
    fig.savefig(OUT.parent / "panels" / f"fig2_{name}.png", dpi=600, bbox_inches=bb)
for a_ in arts:
    a_.set_visible(True)

# ---------------------------------------------------------------- the numbers the figure shows
for k in KINDS:
    print(f"{k:12s} H {R[k]['H']:.4f}" + (f" (mean of seeds {ENS[k]:.4f})" if k in ENS else " = log2 N") +
          f", k {sorted(R[k]['k'].astype(int).tolist(), reverse=True)}, information {np.nanmin(R[k]['I']):.2f} to {np.nanmax(R[k]['I']):.2f}"
          + (f", practice passes it at step {PASS[k]}" if k in PASS else ""))
print("isolated elements:", {k: int((R[k]["k"] == 0).sum()) for k in KINDS})
print(f"practice: mean H {HMEAN[0]:.3f} -> {HMEAN[200]:.3f} (200) -> {HMEAN[500]:.3f} (500) -> {HMEAN[-1]:.3f} ({T_END})")
print(f"beta {BETA:.4f} ms per bit, a {A_MS} ms; slopes A {BETA * CEIL:.2f}, D {BETA * RD['H']:.2f}, t = {T_END} {BETA * HMEAN[-1]:.2f}; "
      f"ratio D/A {RD['H'] / CEIL:.4f}, end/A {HMEAN[-1] / CEIL:.4f}")
print(f"race: theta {THETA}, {UNIT:.2f} ms per unit of 1/v; hub share {RD['p'][hub]:.3f} selected at {T_ARR['hub']:.1f} ms, "
      f"leaf share {RD['p'][weak]:.3f} at {T_ARR['leaf']:.1f} ms, every element of A at {T_ARR['every element']:.1f} ms; "
      f"mean of D {(THETA + np.log(2.0) * RD['H']) * UNIT:.1f} ms")
for bits in (2, 4):
    print(f"T at log2 X = {bits}: {A_MS + BETA * CEIL * bits:.1f} -> {A_MS + BETA * HMEAN[-1] * bits:.1f} ms")
print("written", OUT)
