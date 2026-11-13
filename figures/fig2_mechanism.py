"""Figure 2. The dynamic law as a scheme, T = a + beta H(t) log2 X, from even coupling (A) to hub coupling (D).

The classical laws (Equations 1 to 3) hold the coupling even, H = log2 N, and so describe A alone. The dynamic law covers the
whole path, because practice changes the coupling and with it H(t) and the slope b = beta H(t).

Every structure is read the same way. Node area reads the share p_i (Equation 5, c_i = k_i), link width the coupling
of the pair, and colour log2(p_i N) = log2 N - I(p_i), the bits by which an element's information I(p_i) = -log2 p_i lies below
the even information, on the diverging scale of Park (2026): blue below an even share, pale yellow at 1/N, yellow to red
above it, on +-2.716 bits. Below each structure every element is a block of width p_i and height I(p_i), so the total
area of the blocks is H (Equations 6 and 7), the solid line, with log2 N dashed.

A, even coupling: 4 elements and 4 links in a closed chain (shares 1/4, shortest code 00 01 10 11, H = 2 bits) and the
chain of 20 elements and 20 links (H = log2 20). D, hub coupling, mirrored: the scale-free structure of 20 elements and
the same 4 elements coupled through a hub with a double link (shares 1/2, 1/4, 1/8, 1/8, code 0 10 110 111, H = 1.75
bits). In both four-element structures each codeword is as long as the information of its element, so the mean code length
equals H (Equation 9); the script checks this. Between A and D, unlabelled, the small-world and random structures.
B, H(t) under the practice rules of the full model (code/sim_full.py), whose starting ring is the chain of A: formation
by eta_i k_i, selective elimination by 1/eta per link end from step 200, lambda = 0.5 with g = 1 on three elements; 300 runs,
mean and interquartile band. C, the law at A, at D and at the end of practice (Equation 10).

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
Writes figures/fig2_mechanism.png and prints every number the figure shows.
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
# The palette and the colour scale of Park (2026, Figures 1 to 6), unchanged. One diverging scale for every element:
# blue below an even share, pale yellow at exactly 1/N, yellow to red above it, read as log2(p_i N) in bits on +-2.716.
INK, GREY, GUIDE, RULE, TRACK = "#1A1A1A", "#8A8A8A", "#455A64", "#8C989E", "#E3E7E9"
LINK, EDGE = "#C3CBD0", "#37474F"                                     # links; outlines, stems and arrows
PALE = TRACK
ENS_TONE = {"regular": "#B9C4C9", "small-world": "#93A8B1", "random": "#6E8C99", "scale-free": "#48707F"}
_HOT = plt.get_cmap("hot_r")(np.linspace(0.06, 0.94, 256))
_COOL = plt.get_cmap("Blues")(np.linspace(0.92, 0.06, 256))
CMAP = LinearSegmentedColormap.from_list("hot_div", np.vstack([_COOL, _HOT]))
VLIM = 2.7157337518864986
CNORM = Normalize(-VLIM, VLIM)

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

# ---------------------------------------------------------------- four elements and four links, the coding example
SMALL = {"even": nx.MultiGraph([(0, 1), (1, 2), (2, 3), (3, 0)]),     # a closed chain, every share 1/4
         "hub": nx.MultiGraph([(0, 1), (0, 1), (0, 2), (0, 3)])}     # a hub, with a double link to element 1
CODE = {"even": ("00", "01", "10", "11"), "hub": ("0", "10", "110", "111")}
RS = {k: readings(G) for k, G in SMALL.items()}
for k, c in CODE.items():                    # shortest code: prefix-free, complete, each word as long as its information
    assert not any(a != b and b.startswith(a) for a in c for b in c) and sum(2.0 ** -len(w) for w in c) == 1
    assert all(len(w) == I for w, I in zip(c, RS[k]["I"]))
assert RS["even"]["H"] == 2.0 and RS["hub"]["H"] == 1.75

# ---------------------------------------------------------------- practice from A: the fixed regime of the full model
LAM, T_ELIM, T_END, RUNS = 0.5, 200, 1000, 300
HT_ = np.array([sf.run(N=N, n_start=N, steps=T_END, lam=LAM, t_elim=T_ELIM, seed=s)[0] for s in range(RUNS)])
HMEAN = HT_.mean(0)
PASS = {k: int(np.argmax(HMEAN <= R[k]["H"])) for k in KINDS[1:]}

PHI = max(R[k]["p"].max() for k in KINDS)


def pcol(I, n):
    """Colour of an element of information I among n: log2(p_i n) = log2 n - I(p_i), the bits it sits below the even value."""
    b = np.log2(n) - I
    return CMAP(CNORM(0.0 if abs(b) < 1e-9 else b))


def ink_on(I, n):
    r, g, b, _ = pcol(I, n)
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.45 else INK


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
SPOS = {"even": {0: (-0.6, -0.6), 1: (-0.6, 0.6), 2: (0.6, 0.6), 3: (0.6, -0.6)},
        "hub": {0: (0.0, 0.05), 1: (-0.72, -0.62), 2: (0.0, 0.92), 3: (0.74, -0.62)}}

FW, FH = 6.5, 2.68
Y_TITLE = FH - 0.14


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


def network(ax, kind, size=1.0, sx=1.0):
    """Node area by share p_i, link width by the coupling of the pair, node shade by information I(p_i)."""
    G, Rk = NET[kind], R[kind]
    pos = wide(ax, POS[kind], sx, 1.1)
    for u, v in G.edges():
        ax.plot(*zip(pos[u], pos[v]), color=LINK, lw=link_width(Rk["p"][u], Rk["p"][v], size), solid_capstyle="round", zorder=2)
    for n in range(N):
        if Rk["p"][n] == 0:                                          # no links, no share
            ax.scatter(*pos[n], s=size * 6, color="white", edgecolors=GREY, linewidths=0.45, zorder=5)
            continue
        ax.scatter(*pos[n], s=size * (3 + 560 * Rk["p"][n]), color=[pcol(Rk["I"][n], N)], edgecolors=EDGE,
                   linewidths=0.45, zorder=5)
    return pos


def small_network(ax, kind, size=0.42, sx=1.0):
    """The four-element structures in the same reading, a double link drawn as two strokes, each codeword beside its element."""
    G, Rk = SMALL[kind], RS[kind]
    pos = wide(ax, SPOS[kind], sx, 1.25)
    for (u, v), m in Counter(tuple(sorted(e)) for e in G.edges()).items():
        P, Q = np.array(pos[u]), np.array(pos[v]); d = Q - P; nrm = np.array([-d[1], d[0]]) / np.hypot(*d)
        for off in (np.linspace(-0.10, 0.10, m) if m > 1 else [0.0]):
            ax.plot(*zip(P + off * nrm, Q + off * nrm), color=LINK, lw=link_width(Rk["p"][u], Rk["p"][v], size),
                    solid_capstyle="round", zorder=2)
    cen = np.mean([pos[n] for n in pos], axis=0)
    for n in sorted(G.nodes()):
        ax.scatter(*pos[n], s=size * (3 + 560 * Rk["p"][n]), color=[pcol(Rk["I"][n], 4)], edgecolors=EDGE, linewidths=0.45, zorder=5)
        out = np.array(pos[n]) - cen; out = out / (np.hypot(*out) or 1)
        if kind == "hub" and n == 0:
            out = np.array([-1.0, 0.25]) / np.hypot(1.0, 0.25)
        r_ = 0.36 + 0.5 * np.sqrt(Rk["p"][n])
        ax.text(pos[n][0] + r_ * out[0], pos[n][1] + r_ * out[1], CODE[kind][n], fontsize=6.6, ha="center", va="center", zorder=6)
    return pos


YMAX = 6.0


def info_blocks(ax, Rk, n, mirror=False, order=None, codes=None, ylabel=True):
    """Blocks of width p_i and height I(p_i), in order of rising share from the outer edge inward (ties in the
    given order). The right-hand panels are mirrored, so their shares run from right to left."""
    cum, top = 0.0, {}
    for i in (np.argsort(Rk["p"], kind="stable") if order is None else order):
        ax.add_patch(Rectangle((cum, 0), Rk["p"][i], Rk["I"][i], facecolor=pcol(Rk["I"][i], n), edgecolor=EDGE, lw=0.25, zorder=3))
        if codes is not None:                                       # a narrow block holds its codeword upright, below the H line
            up = Rk["p"][i] < 0.2; y_ = Rk["I"][i] / 2
            if up and Rk["I"][i] > Rk["H"]:
                y_ = min(y_, (Rk["H"] - 0.12) / 2)
            ax.text(cum + Rk["p"][i] / 2, y_, codes[i], rotation=90 if up else 0, fontsize=6.4, ha="center",
                    va="center", color=ink_on(Rk["I"][i], n), zorder=7)
        top[i] = (cum + Rk["p"][i] / 2, Rk["I"][i]); cum += Rk["p"][i]
    ax.plot([0, 1], [np.log2(n)] * 2, color=RULE, lw=0.7, ls=(0, (2.4, 1.6)), zorder=6)
    ax.plot([0, 1], [Rk["H"]] * 2, color=INK, lw=1.3, zorder=8)
    ax.set_xlim((1, 0) if mirror else (0, 1)); ax.set_ylim(0, YMAX); ax.set_yticks([0, 2, 4])
    if n == N:
        ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    else:
        ax.set_xticks([0, 1]); ax.set_xticklabels(["0", "1"])
    ax.spines["top"].set_visible(False)
    ax.spines["left" if mirror else "right"].set_visible(False)
    if mirror:
        ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
    if ylabel:
        ax.set_ylabel(r"information $I_i$, bits", labelpad=1.5)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(rf"shares, $N$ = {n}", labelpad=0.8)
    return top


def guide(axA, xyA, axB, xyB, rad=0.0):
    fig.add_artist(ConnectionPatch(xyA=xyA, coordsA=axA.transData, xyB=xyB, coordsB=axB.transData, arrowstyle="-|>",
                                   mutation_scale=6, lw=0.7, color=GUIDE, ls=(0, (2.0, 1.6)), zorder=9,
                                   shrinkA=4, shrinkB=1, connectionstyle=f"arc3,rad={rad}"))


def fig_arrow(x0, x1, y):
    fig.add_artist(FancyArrowPatch((x0 / FW, y / FH), (x1 / FW, y / FH), transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=8, lw=0.9, color=EDGE, ls=(0, (2.4, 1.8)), zorder=9))


def top_bottom_box(ax):
    for sp in ("left", "right"):
        ax.spines[sp].set_visible(False)


def scope(x0, x1, y, label, tick=0.07):
    """A bar under the panels it covers, open upward at both ends, its label set in the bar."""
    fig.add_artist(Line2D(np.array([x0, x0, x1, x1]) / FW, np.array([y + tick, y, y, y + tick]) / FH,
                          transform=fig.transFigure, color=INK, lw=0.7))
    fig.text((x0 + x1) / 2 / FW, y / FH, label, fontsize=7.2, ha="center", va="center",
             bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none"))


def ftext(x, y, s, **kw):
    return fig.text(x / FW, y / FH, s, **kw)


KEY = dict(frameon=False, fontsize=6.4, handlelength=1.5, handletextpad=0.4, borderaxespad=0.1, labelspacing=0.2)
SOLID = Line2D([], [], color=INK, lw=1.3); DASHED = Line2D([], [], color=GREY, lw=0.7, ls=(0, (2.4, 1.6)))

# ---------------------------------------------------------------- columns
X_A4, W_4, W_20, GAP = 0.32, 0.70, 0.94, 0.10
X_A20 = X_A4 + W_4 + GAP
X_D4 = FW - X_A4 - W_4; X_D20 = X_D4 - GAP - W_20
XC_A = (X_A4 + X_A20 + W_20) / 2; XC_D = (X_D20 + X_D4 + W_4) / 2
Y_P, H_P, H_N = 0.92, 0.56, 0.84
Y_N = Y_P + H_P + 0.12
SX = 1.2                                                              # sideways stretch of every network drawing

fig = plt.figure(figsize=(FW, FH))


def code_in(ax, lines, x0=0.03, y0=0.97, dy=0.18, fs=6.4):
    """The sum for H written in the empty top of an information panel, continuation lines set under the first term.
    Drawn on the figure above the guide arrows, so the arrows pass behind the lines."""
    rend = fig.canvas.get_renderer(); bb = ax.get_position()
    t_ = fig.text(0, 0, r"$H\ $", fontsize=fs); ind = t_.get_window_extent(rend).width / fig.dpi / (bb.width * FW); t_.remove()
    texts = [fig.text(x0 + (ind if k else 0.0), y0 - k * dy, s_, transform=ax.transAxes, fontsize=fs, ha="left",
                      va="center", zorder=12) for k, s_ in enumerate(lines)]
    for t in texts:                                                   # a white ground under each line, all grounds below all lines
        (x_0, y_0), (x_1, y_1) = fig.transFigure.inverted().transform(t.get_window_extent(rend).padded(0.8 * fig.dpi / 72).get_points())
        fig.add_artist(Rectangle((x_0, y_0), x_1 - x_0, y_1 - y_0, transform=fig.transFigure, fc="white", ec="none", zorder=11))


# ---------------------------------------------------------------- A, the even structures, 4 and 20 elements
heading("A", "Novice, even coupling", 0.04, "left")
p4 = ax_in(X_A4, Y_P, W_4, H_P); top4 = info_blocks(p4, RS["even"], 4, order=[0, 1, 2, 3], codes=CODE["even"])
s4 = ax_in(X_A4 + W_4 / 2 - 0.37, Y_N + 0.10, 0.74, 0.64); ps4 = small_network(s4, "even", sx=SX)
for n in (0, 3):                                                      # the two lower elements, each to its block
    guide(s4, ps4[n], p4, top4[n])
code_in(p4, [r"$H = \frac{1}{4}· 2+\frac{1}{4}· 2$", r"$+\,\frac{1}{4}· 2+\frac{1}{4}· 2$", r"$=\,2$ bits"], x0=0.07)
RA = R["regular"]
by_x = sorted(range(N), key=lambda n: LOOP[n][0])                    # equal shares, ties ordered by position
pA = ax_in(X_A20, Y_P, W_20, H_P); topA = info_blocks(pA, RA, N, order=by_x, ylabel=False)
pA.text(0.5, CEIL + 0.3, rf"$H$ = {CEIL:.3f}", fontsize=6.6, ha="center", va="bottom")
aA = ax_in(X_A20 + W_20 / 2 - 0.54, Y_N, 1.00, H_N); pa = network(aA, "regular", sx=SX)
low = sorted((n for n in range(N) if LOOP[n][1] < 0), key=lambda n: LOOP[n][0])   # two lower elements, near the ends of the row
for n in (low[1], low[-1]):
    guide(aA, pa[n], pA, topA[n])

# ---------------------------------------------------------------- D, the hub structures, 20 and 4 elements, mirrored
heading("D", "Expert, hub coupling", FW - 0.04, "right")
RD = R["scale-free"]; pos_d = POS["scale-free"]
leaves = np.flatnonzero(RD["k"] == RD["k"].min())
weak = int(max(leaves, key=lambda n: pos_d[n][0] - pos_d[n][1]))     # a leaf toward the lower right
pD = ax_in(X_D20, Y_P, W_20, H_P)
topD = info_blocks(pD, RD, N, mirror=True, order=sorted(range(N), key=lambda n: (RD["p"][n], n != weak)), ylabel=False)
pD.text(0.60, CEIL + 0.3, rf"$H$ = {RD['H']:.3f}", fontsize=6.6, ha="center", va="bottom")
aD = ax_in(X_D20 + W_20 / 2 - 0.42, Y_N, 1.00, H_N); pd_ = network(aD, "scale-free", sx=SX)
hub = int(np.argmax(RD["p"]))
guide(aD, pd_[hub], pD, topD[hub], rad=0.35)
guide(aD, pd_[weak], pD, topD[weak])
q4 = ax_in(X_D4, Y_P, W_4, H_P); topq = info_blocks(q4, RS["hub"], 4, mirror=True, order=[3, 2, 1, 0], codes=CODE["hub"])
t4 = ax_in(X_D4 + W_4 / 2 - 0.37, Y_N + 0.10, 0.74, 0.64); pt4 = small_network(t4, "hub", sx=SX)
guide(t4, pt4[0], q4, topq[0], rad=0.25)                              # the hub and its weakest element, each to its block
guide(t4, pt4[3], q4, topq[3])
code_in(q4, [r"$H = \frac{1}{2}· 1+\frac{1}{4}· 2$", r"$+\,\frac{1}{8}· 3+\frac{1}{8}· 3$", r"$=\,1.75$ bits"])

# ---------------------------------------------------------------- B, H(t) from A, with the levels of A and D
X_B, W_B, X_C, W_C, Y_T, H_T = 2.40, 0.54, 3.16, 0.94, Y_N + 0.14, 0.72
XM = (X_A20 + W_20 + X_D20) / 2                                       # the middle, between A and D
head("B", r"Practice, $H(t)$", X_B + W_B / 2)
bB = ax_in(X_B, Y_T, W_B, H_T)
t = np.arange(T_END + 1)
bB.fill_between(t, np.percentile(HT_, 25, axis=0), np.percentile(HT_, 75, axis=0), color=PALE, lw=0, zorder=1)
bB.plot(t, HMEAN, color=INK, lw=1.4, zorder=4)
for k, lab in (("regular", "A"), ("scale-free", "D")):
    bB.axhline(R[k]["H"], color=RULE, lw=0.6, ls=(0, (1.0, 1.4)), zorder=2)
    bB.text(T_END * 1.03, R[k]["H"], lab, fontsize=7, fontweight=WEIGHT, ha="left", va="center", clip_on=False)
for k in ("small-world", "random"):                                  # where practice passes the unlabelled structures
    bB.plot([PASS[k]], [R[k]["H"]], "o", ms=2.8, color=ENS_TONE[k], mec=EDGE, mew=0.5, zorder=6)
bB.plot([PASS["scale-free"]], [RD["H"]], "o", ms=3.2, color=ENS_TONE["scale-free"], mec=EDGE, mew=0.5, zorder=6)
bB.set_xlim(0, T_END); bB.set_ylim(2.85, 4.48); bB.set_yticks([3.0, 3.5, 4.0]); bB.set_xticks([0, 500, 1000])
top_bottom_box(bB)
bB.set_xlabel("practice step $t$", labelpad=0.8); bB.set_ylabel(r"$H(t)$, bits", labelpad=1.5)

# ---------------------------------------------------------------- C, the law at A, at D and at the end of practice
head("C", r"Slope, $b = \beta H(t)$", X_C + W_C / 2)
cC = ax_in(X_C, Y_T, W_C, H_T)
xb = np.array([0.0, 4.0])
rows = [(CEIL, dict(color=ENS_TONE["small-world"], lw=1.3), "A"), (RD["H"], dict(color=ENS_TONE["scale-free"], lw=1.3), "D"),
        (HMEAN[-1], dict(color=INK, lw=1.6), f"$t$ = {T_END}")]
for H_, kw, lab in rows:
    cC.plot(xb, A_MS + BETA * H_ * xb, label=f"{lab}, {BETA * H_:.1f}", **kw)
cC.annotate("", xy=(4.0, A_MS + BETA * HMEAN[-1] * 4.0), xytext=(4.0, A_MS + BETA * CEIL * 4.0),
            arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=0.7, color=INK), annotation_clip=False)
cC.legend(loc="upper left", fontsize=6.2, frameon=False, handlelength=0.9, handletextpad=0.3, borderaxespad=0.02,
          labelspacing=0.1)                                          # slope b of each line, in ms per bit
cC.text(-0.08, A_MS, "$a$", fontsize=7, ha="right", va="center", clip_on=False)
cC.set_xlim(0, 4.0); cC.set_xticks([0, 2, 4]); cC.set_ylim(0, 445); cC.set_yticks([0, 200, 400])
top_bottom_box(cC); cC.yaxis.tick_right(); cC.yaxis.set_label_position("right")
cC.set_xlabel(r"$\log_2 X$, bits", labelpad=0.8); cC.set_ylabel(r"$T$, ms", labelpad=1.5)

# ---------------------------------------------------------------- between A and D, small, the small-world and random structures
Y_M, S_M = Y_P + H_P / 2 - 0.22, 0.44
sz = (S_M / H_N) ** 2
aS = ax_in(XM - 0.55 - S_M / 2, Y_M, S_M, S_M); network(aS, "small-world", size=sz)
aR = ax_in(XM + 0.55 - S_M / 2, Y_M, S_M, S_M); network(aR, "random", size=sz)
y_mid = Y_M + S_M / 2
ftext(XM, y_mid, r"$\cdots$", fontsize=12, ha="center", va="center")
fig_arrow(X_A20 + W_20 + 0.06, XM - 0.55 - S_M / 2 - 0.02, y_mid); fig_arrow(XM + 0.55 + S_M / 2 + 0.02, X_D20 - 0.06, y_mid)
ftext(XM, y_mid - 0.09, "practice $t$", fontsize=6.8, ha="center", va="top")

# ---------------------------------------------------------------- below A the classical laws, below D the slope of the dynamic law
EQ = dict(fontsize=7.0, ha="left", va="center")
y2, y3 = 0.555, 0.41
ftext(0.10, (y2 + y3) / 2, r"$RT = a + b\,\log_2 n$", **EQ)
ftext(1.02, y2, r"$MT = a + b\,\log_2(2D/W)$", **EQ); ftext(1.02, y3, r"$W_e = a + b\,(D/MT)$", **EQ)
ftext(XC_D, y2, r"$b = \beta H(t) < \beta\,\log_2 N$", fontsize=7.4, ha="center", va="center")
cb = ax_in(XC_D - 0.62, y3 - 0.025, 0.86, 0.05)                         # the colour scale of every node and block
cbar = mpl.colorbar.ColorbarBase(cb, cmap=CMAP, norm=CNORM, orientation="horizontal", ticks=[-2, -1, 0, 1, 2])
cb.tick_params(labelsize=6.0, length=1.6, width=0.5, pad=1.0); cbar.outline.set_linewidth(0.5)
cb.text(1.04, 0.5, r"$\log_2(p_i N)$, bits", transform=cb.transAxes, fontsize=6.4, ha="left", va="center")

# ---------------------------------------------------------------- the law in the middle, its terms named
y_eq = 0.62
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

# ---------------------------------------------------------------- what each law covers: the classical laws A alone, the dynamic law A to D
scope(0.06, X_A20 + W_20 + 0.04, 0.26, "classical laws, static")
scope(0.06, FW - 0.06, 0.10, r"dynamic law, $H(t)$ from A to D")

# the trend lines of B and C set narrow, as in figures/style.py (bolden): lines of at least 0.95 pt in axes with ticks
FACTOR = 1.4 if BOLD else 0.75
for ax_ in fig.axes:
    if len(ax_.get_xticks()) or len(ax_.get_yticks()):
        for ln in ax_.get_lines():
            if ln.get_linestyle() not in ("None", "", " ") and ln.get_linewidth() >= 0.95:
                ln.set_linewidth(ln.get_linewidth() * FACTOR)
        leg = ax_.get_legend()
        for h_ in (leg.legend_handles if leg is not None else []):
            if hasattr(h_, "set_linewidth") and h_.get_linestyle() not in ("None", "", " ") and h_.get_linewidth() >= 0.95:
                h_.set_linewidth(h_.get_linewidth() * FACTOR)
fig.savefig(ROOT / "figures" / "fig2_mechanism.png", dpi=600)

# each lettered panel on its own: the artists of the figure are sorted by where their centres lie, the others are hidden
# while the panel is saved, and the crop is the union of what remains. The centre strip, the law and practice between A
# and D, has no letter and is written as fig2_centre; the bar of the dynamic law spans the figure and is left out.
from matplotlib.transforms import Bbox  # noqa: E402
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.dpi_scale_trans.inverted()
X_END_A, X_BEG_D, Y_LOW = X_A20 + W_20 + 0.10, X_D20 - 0.04, Y_N - 0.32


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


arts = [a_ for a_ in list(fig.axes) + list(fig.texts) + list(fig.artists) + list(fig.lines) + list(fig.patches)
        if a_.get_visible()]
boxes = [a_.get_window_extent(rend).transformed(inv) if not hasattr(a_, "get_tightbbox") or a_.get_tightbbox(rend) is None
         else a_.get_tightbbox(rend).transformed(inv) for a_ in arts]
where = [panel_of(bb) for bb in boxes]
for name in ("A", "B", "C", "D", "centre"):
    keep = [i for i, w in enumerate(where) if w == name]
    for i, a_ in enumerate(arts):
        a_.set_visible(i in keep)
    bb = Bbox.union([boxes[i] for i in keep]).padded(0.04)
    fig.savefig(ROOT / "figures" / "panels" / f"fig2_{name}.png", dpi=600, bbox_inches=bb)
for a_ in arts:
    a_.set_visible(True)

for k in KINDS:
    print(f"{k:12s} H {R[k]['H']:.4f}" + (f" (mean of seeds {ENS[k]:.4f})" if k in ENS else " = log2 N") +
          f", k {sorted(R[k]['k'].astype(int).tolist(), reverse=True)}, information {np.nanmin(R[k]['I']):.2f} to {np.nanmax(R[k]['I']):.2f}"
          + (f", practice passes it at step {PASS[k]}" if k in PASS else ""))
for k in SMALL:
    print(f"4 elements, {k:4s} k {RS[k]['k'].astype(int).tolist()}, shares {RS[k]['p'].tolist()}, codes {CODE[k]}, H {RS[k]['H']:.3f}")
print("isolated elements:", {k: int((R[k]["k"] == 0).sum()) for k in KINDS})
print(f"practice: mean H {HMEAN[0]:.3f} -> {HMEAN[200]:.3f} (200) -> {HMEAN[500]:.3f} (500) -> {HMEAN[-1]:.3f} ({T_END})")
print(f"beta {BETA:.4f} ms per bit, a {A_MS} ms; slopes A {BETA * CEIL:.2f}, D {BETA * RD['H']:.2f}, t = {T_END} {BETA * HMEAN[-1]:.2f}; "
      f"ratio D/A {RD['H'] / CEIL:.4f}, end/A {HMEAN[-1] / CEIL:.4f}")
for bits in (2, 4):
    print(f"T at log2 X = {bits}: {A_MS + BETA * CEIL * bits:.1f} -> {A_MS + BETA * HMEAN[-1] * bits:.1f} ms")
