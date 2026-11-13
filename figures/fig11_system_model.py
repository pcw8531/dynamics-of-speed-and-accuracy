"""Figure 11 of the manuscript: a systems model of the dynamic law (A) and the law through practice (B), the mechanism
of Figures 2 to 10 in one diagram.

A, the law as a causal-loop diagram. On the left, the initial state: N elements share coupling evenly, k_i(0) the same
for every element, so H(0) = log2 N, drawn as a small ring, and practice directed at some of them, g_i, enters as an
input to a targeted element and sets its fitness. The coupling held by each element, k_i(t), is the state of the
system. Formation and recruitment add to it and elimination removes from it. Arrows carry influence, each with the sign
of its effect. Coupling feeds formation through preferential attachment, the reinforcing loop R, and elimination through
the links it holds, the balancing loop B, and fitness, eta_i = 1 + lambda g_i, weights formation toward fit elements and
elimination away from them. Elimination lowers H only because it is selective, by the weight 1/eta_i. The shares of
coupling give the averaged information H(t), which sets the slope of the time of choosing and aiming, and with the
correlations among joined elements they give the relative spread of force, which sets the slope of the
impulse-variability equation. A dotted frame holds the law, and both converge on it: H(t) sets its slope, and the
impulse-variability equation gives the spread the initial adjustment leaves, whose error, E_0, sets the bits
X = 2 E_0 / W that the corrections remove at beta H(t) each (Equation 20). Solid arrows are the rules the simulations
run, dashed arrows links derived in the text and not simulated (Equations 20 and 22). Grey labels mark what earlier
accounts take as given: the classical laws read the law at one step, and accounts of coordination describe the shares
and correlations at one level of skill. No data enter A.
B, left, at one step: the law T = a + beta H(t) log2 X at successive steps of practice, with H(t) the mean of the
fixed condition of the full model (sim_out/sim_full_results.npz), a = 12.8 ms and beta = 94.7 / log2 20 ms per command
bit as in Figures 2, 4 and 5. Each step is a straight line, a classical law, and the line at t = 0 is the case of the
classical laws; each later line takes the colour that log2 N - H(t) bits has on the scale of Figure 2. B, right, over
practice: the same law across practice for tasks of 1, 4 and 7 bits, with the time saved by step 3,000 at 7 bits and the
time axis repeated on the right. One network sits in each plot, with a dotted arrow to its point, as in Figures 4B and
5B: in the left plot the even ring of step 0, on the line of step 0, and in the right plot the coupling of the simulated
performer of Figures 6 to 9 at step 1,000 (the run of code/sim_static_cases.py, seed 28, H = 3.006 bits, rebuilt here
with code/sim_full.py and checked against the stored shares), on the 7-bit course. Both are drawn in the reading of
Figure 2, node area by share, node colour by log2(p_i N) and link width and shade by the number of links between the
pair: the ring as a closed chain of another shape than in Figures 2A, 4B and 5B, and the performer's network with its
three hubs held apart at the centre and the other elements placed by a spring layout (networkx), its three elements
without links as small open circles below.
A and B stand side by side. Palette and type from figures/style.py. Each panel is also written on its own to figures/panels.
Run from the repository root: python figures/fig11_system_model.py sim_out
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Arc, Circle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, EDGE, RULE, GREY, TONE, share_colour, tb_axes, save_panels, inset_at, loop_layout  # noqa: E402
sys.path.insert(0, str(ROOT / "code"))
import sim_full as sf  # noqa: E402

src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "sim_out"
d = np.load(src / "sim_full_results.npz")
SC = np.load(src / "sim_static_cases.npz")                            # the simulated performer of Figures 6 to 9
H = d["fixed"]; N, A_MS = 20, 12.8
BETA = 94.7 / np.log2(N)
CEIL = np.log2(N)
INFO = TONE["scale-free"]                                            # the colour of influence arrows

FW, FH = 6.5, 2.24
S = Sheet(FW, FH)
XL0, XR0, Y0, Y1 = -0.52, 3.60, 0.19, 2.17                            # panel A, drawn in inches
SCALE = 0.870                                                        # and set smaller than drawn
F = 0.93                                                             # type of A, a little smaller to match
AX0, AY = 0.02, 0.17
AWID, AH = (XR0 - XL0) * SCALE, (Y1 - Y0) * SCALE
ax = S.ax(AX0, AY, AWID, AH)
ax.set_xlim(XL0, XR0); ax.set_ylim(Y0, Y1); ax.axis("off")


# ---------------------------------------------------------------- drawing primitives of the diagram
def node(xc, yc, w, h, lines, kind="process"):
    """A variable of the system: the state (coupling) in a filled box, a process in a light rounded box, a line of the
    law in a box with a heavy border, a plain variable as text alone."""
    if kind == "state":
        ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.03",
                                    fc=share_colour(0.0), ec=INK, lw=1.2, zorder=4))
    elif kind == "process":
        ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.08",
                                    fc="white", ec=EDGE, lw=0.7, zorder=4))
    elif kind == "law":
        ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.04",
                                    fc="white", ec=INK, lw=1.2, zorder=4))
    n = len(lines); dy = 0.13
    for j, ln in enumerate(lines):
        grey = kind == "law" and j > 0
        ax.text(xc, yc + (n - 1) * dy / 2 - j * dy, ln, fontsize=F * (5.9 if grey else (6.6 if j == 0 else 6.4)),
                ha="center", va="center", zorder=5, color=EDGE if grey else INK, style="italic" if grey else "normal")


def info(p0, p1, sign=None, rad=0.0, at=0.86, off=(0.05, 0.0), dashed=False):
    """An arrow of influence, with the sign of its effect near the head."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=6.5, lw=0.8, color=INFO, zorder=6,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=1.5, shrinkB=1.5,
                                 linestyle=(0, (2.4, 1.6)) if dashed else "-"))
    if sign:
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float); mid = (p0 + p1) / 2; dv = p1 - p0
        ctrl = mid + rad * np.array([dv[1], -dv[0]])
        q = (1 - at) ** 2 * p0 + 2 * (1 - at) * at * ctrl + at ** 2 * p1
        ax.text(q[0] + off[0], q[1] + off[1], sign, fontsize=F * 7.6, color=INFO, ha="center", va="center", zorder=7)


def loop(x, y, letter, ccw=True):
    """The mark of a feedback loop: a circular arrow around R (reinforcing) or B (balancing)."""
    ax.add_patch(Arc((x, y), 0.19, 0.19, theta1=60, theta2=340, color=EDGE, lw=0.7, zorder=6))
    th = np.deg2rad(340 if ccw else 60); tip = np.array([x + 0.095 * np.cos(th), y + 0.095 * np.sin(th)])
    tan = np.array([-np.sin(th), np.cos(th)]) * (1 if ccw else -1); rad_ = np.array([np.cos(th), np.sin(th)])
    ax.add_patch(Polygon([tip + 0.034 * tan, tip - 0.02 * tan + 0.022 * rad_, tip - 0.02 * tan - 0.022 * rad_],
                         closed=True, fc=EDGE, ec=EDGE, lw=0.3, zorder=6))
    ax.text(x, y, letter, fontsize=F * 6.8, ha="center", va="center", zorder=7)


def note(x, y, txt, ha="center"):
    ax.text(x, y, txt, fontsize=F * 5.9, ha=ha, va="center", color=GREY, style="italic", linespacing=1.05, zorder=5)


def label(x, y, txt, fs, **kw):
    ax.text(x, y, txt, fontsize=F * fs, **kw)


# ---------------------------------------------------------------- A, the loops of practice and what the shares set
CX, CY = 1.42, 1.18                                                  # the state of the system
node(CX, CY, 0.78, 0.34, ["coupling", r"$k_i(t)$"], "state")
FXY, EXY = (0.80, 1.94), (0.80, 0.42)                                # formation above, elimination below
node(*FXY, 0.72, 0.30, ["formation", r"$\Pi_i \propto \eta_i\,k_i$"])
node(*EXY, 0.72, 0.30, ["elimination", r"weight $1/\eta_i$"])
RXY = (1.78, 2.02)
node(*RXY, 0.92, 0.22, [r"recruitment, $m$ links"])

# the initial state above fitness: even coupling over the elements, and practice directed at some of them as an input
IX, IY, IR = 0.0, 1.56, 0.17
label(IX, 2.10, "initial state", 6.6, ha="center", va="center")
label(IX, 1.97, r"$H(0) = \log_2 N$", 6.3, ha="center", va="center")
ang = np.deg2rad(90 + 45 * np.arange(8))
pts = np.c_[IX + IR * np.cos(ang), IY + IR * np.sin(ang)]
for j in range(8):
    q0, q1 = pts[j], pts[(j + 1) % 8]
    ax.plot([q0[0], q1[0]], [q0[1], q1[1]], color=RULE, lw=0.6, zorder=3)
for j, (xq, yq) in enumerate(pts):
    ax.add_patch(Circle((xq, yq), 0.034, fc=share_colour(0.0), ec=EDGE, lw=0.5, zorder=4))
TGT = pts[2]                                                         # the targeted element, on the left of the ring
ax.add_patch(Circle(TGT, 0.060, fc="none", ec=INK, lw=0.7, zorder=5))
UX = TGT[0] - 0.27
ax.add_patch(Circle((UX, TGT[1]), 0.026, fc=INK, ec=INK, zorder=5))
zz = np.linspace(0, 1, 60)                                           # the input, drawn as a wavy arrow as in control
ax.plot(UX + 0.035 + (0.17 - 0.035) * zz, TGT[1] + 0.025 * np.sin(2 * np.pi * 3 * zz), color=INK, lw=0.6, zorder=5)
ax.add_patch(FancyArrowPatch((UX + 0.155, TGT[1]), (TGT[0] - 0.066, TGT[1]), arrowstyle="-|>", mutation_scale=5.5,
                             lw=0.6, color=INK, shrinkA=0, shrinkB=0, zorder=5))
label(UX, TGT[1] + 0.07, r"$g_i$", 6.8, ha="center", va="bottom")
FT = 1.08                                                            # fitness, set by the input
label(IX, FT, "fitness", 6.6, ha="center", va="center")
label(IX, FT - 0.13, r"$\eta_i = 1 + \lambda g_i$", 6.3, ha="center", va="center")
note(IX, FT - 0.32, "quality of\npractice")
info((IX, IY - IR - 0.05), (IX, FT + 0.08), "+", at=0.55, off=(0.06, 0.0))
info((IX + 0.22, FT + 0.08), (FXY[0] - 0.30, FXY[1] - 0.16), "+", at=0.80, off=(-0.07, 0.02))   # fitness weights both flows
info((IX + 0.22, FT - 0.40), (EXY[0] - 0.30, EXY[1] + 0.16), "−", at=0.78, off=(-0.07, -0.02))

# the reinforcing loop of preferential attachment: formation adds coupling, coupling draws formation
F1, C1 = (FXY[0] + 0.22, FXY[1] - 0.16), (CX - 0.30, CY + 0.18)
C2, F2 = (CX - 0.06, CY + 0.18), (FXY[0] + 0.41, FXY[1] - 0.02)
info(F1, C1, "+", rad=0.55, at=0.80, off=(-0.09, 0.00))
info(C2, F2, "+", rad=0.55, at=0.82, off=(0.08, 0.05))
loop(np.mean([F1[0], C1[0], C2[0], F2[0]]) + 0.02, np.mean([F1[1], C1[1], C2[1], F2[1]]), "R")
# the balancing loop of elimination: coupling draws elimination, elimination removes coupling
C3, E3 = (CX - 0.30, CY - 0.18), (EXY[0] + 0.22, EXY[1] + 0.16)
E4, C4 = (EXY[0] + 0.41, EXY[1] + 0.02), (CX - 0.06, CY - 0.18)
info(C3, E3, "+", rad=0.55, at=0.80, off=(-0.09, 0.00))
info(E4, C4, "−", rad=0.55, at=0.84, off=(0.09, -0.04))
loop(np.mean([C3[0], E3[0], E4[0], C4[0]]) + 0.02, np.mean([C3[1], E3[1], E4[1], C4[1]]), "B", ccw=False)
info((RXY[0] - 0.10, RXY[1] - 0.12), (CX + 0.16, CY + 0.18), "+", at=0.76, off=(0.08, 0.02))

# the shares and correlations, and the one law they set, its two lines joined by the error the initial adjustment leaves
SXY, RJ = (2.26, 1.18), (2.12, 0.42)
label(*SXY, r"shares $p_i$", 6.6, ha="center", va="center")
label(*RJ, r"correlations $r_{ij}$", 6.4, ha="center", va="center")
info((CX + 0.41, CY), (SXY[0] - 0.29, SXY[1]))
info((CX + 0.30, CY - 0.18), (RJ[0] - 0.22, RJ[1] + 0.08), "+", rad=0.15, dashed=True, at=0.80, off=(-0.07, -0.04))
note(2.12, 0.80, "coordination,\none level of skill")
LX, FY0, FY1 = 3.10, 1.325, 1.655                                     # the column of the law; a dotted frame holds its key equation
ax.add_patch(FancyBboxPatch((LX - 0.48, FY0), 0.96, FY1 - FY0, boxstyle="round,pad=0,rounding_size=0.05",
                            fc="white", ec=INK, lw=0.8, ls=(0, (1.0, 1.6)), zorder=4))
label(LX, 2.075, r"$H(t) = -\Sigma_i\,p_i \log_2 p_i$", 6.4, ha="center", va="center", zorder=5)
label(LX, 1.945, "averaged information", 5.9, ha="center", va="center", color=EDGE, style="italic", zorder=5)
label(LX, 1.545, r"$T = a + \beta\,H(t)\,\log_2 X$", 7.4, ha="center", va="center", zorder=5)
label(LX, 1.405, "choosing and aiming", 5.9, ha="center", va="center", color=EDGE, style="italic", zorder=5)
label(LX - 0.05, 0.865, r"$W_e = a + b(t)\,D/MT$", 6.6, ha="center", va="center", zorder=5)
label(LX - 0.05, 0.725, "initial adjustment", 5.9, ha="center", va="center", color=EDGE, style="italic", zorder=5)
HA, XA = LX + 0.03, LX + 0.37                                        # above the H(t) and under the X of the first line
info((HA, 1.885), (HA, FY1), "+", at=0.55, off=(0.07, 0.0))
info((XA, 0.94), (XA, FY0), "+", at=0.80, off=(0.065, 0.0), dashed=True)
label(XA - 0.06, 1.205, "error left,", 5.9, ha="right", va="center", color=EDGE, style="italic", zorder=5)
label(XA - 0.06, 1.075, r"$X = 2E_0/W$", 6.2, ha="right", va="center", zorder=5)
FSXY = (2.70, 0.34)                                                  # the spread of force, below the force line
label(FSXY[0], FSXY[1] + 0.02, r"$\sigma_F / F$", 6.4, ha="center", va="bottom")
label(FSXY[0], FSXY[1] - 0.01, "spread of force", 5.9, ha="center", va="top", color=EDGE)
info((SXY[0] + 0.12, SXY[1] + 0.11), (LX - 0.40, 2.04), rad=-0.15)
info((SXY[0] + 0.12, SXY[1] - 0.11), (FSXY[0] - 0.06, FSXY[1] + 0.14), rad=0.10, dashed=True)  # Equation 22, not simulated
info((RJ[0] + 0.32, RJ[1] - 0.02), (FSXY[0] - 0.13, FSXY[1] + 0.06), "+", dashed=True, at=0.55, off=(0.0, 0.07))
info((FSXY[0] + 0.06, FSXY[1] + 0.13), (FSXY[0] + 0.06, 0.665), "+", at=0.50, off=(0.07, 0.0), dashed=True)
note(LX + 0.48, 0.30, "classical laws,\none step", ha="right")

# ---------------------------------------------------------------- B, left: the law at successive steps of practice
BX, BX2, BY, BW, BW2, BH = 3.98, 5.18, 0.36, 1.04, 0.94, 1.52     # B, left at one step and right over practice
axB = S.ax(BX, BY, BW, BH)
x = np.linspace(0, 7, 50)
STEPS = (3000, 2000, 1500, 1000, 700, 400, 200, 100, 0)
for tt in STEPS:
    z_ = A_MS + BETA * H[tt] * x
    c_ = INK if tt == 0 else share_colour(CEIL - H[tt])
    axB.plot(x, z_, color=c_, lw=0.9 if tt == 0 else 0.7, zorder=3 if tt == 0 else 2, solid_capstyle="round")
axB.text(6.9, A_MS + BETA * H[0] * 7 + 12, "$t$ = 0", fontsize=5.8, ha="right", va="bottom")
axB.text(6.95, 120, "$t$ = 3,000", fontsize=5.8, ha="right", va="top", color=EDGE)
axB.text(0.32, -16, "$a$", fontsize=6.6, ha="left", va="center")
tb_axes(axB, xticks=[0, 3.5, 7], xlim=(0, 7), yticks=[0, 350, 700], ylim=(-45, 1000))
axB.set_xlabel(r"$\log_2 X$, bits", labelpad=1); axB.set_ylabel("$T$, ms", labelpad=1.5)

# ---------------------------------------------------------------- B, right: the course of practice on the same time axis
axC = S.ax(BX2, BY, BW2, BH)
t = np.arange(len(H))
STY = {7: dict(color=INK, ls="-"), 4: dict(color=TONE["scale-free"], ls=(0, (4.0, 2.0))),
       1: dict(color=TONE["small-world"], ls=(0, (1.2, 1.6)))}
for bits, st in STY.items():
    axC.plot(t, A_MS + BETA * H * bits, lw=0.9, **st)
axC.text(120, A_MS + BETA * H[0] * 7 + 16, "7 bits", fontsize=5.8, ha="left", va="bottom")
axC.text(120, A_MS + BETA * H[0] * 4 + 16, "4 bits", fontsize=5.8, ha="left", va="bottom", color=TONE["scale-free"])
axC.text(120, A_MS + BETA * H[0] * 1 + 16, "1 bit", fontsize=5.8, ha="left", va="bottom", color=TONE["small-world"])
save = {b_: BETA * (H[0] - H[3000]) * b_ for b_ in (1, 7)}
axC.annotate("", xy=(2890, A_MS + BETA * H[3000] * 7), xytext=(2890, A_MS + BETA * H[0] * 7),
             arrowprops=dict(arrowstyle="-|>", mutation_scale=5, lw=0.6, color=EDGE))
axC.text(2820, 420, f"{save[7]:.0f} ms\nsaved", fontsize=5.8, ha="right", va="center", linespacing=1.0)
tb_axes(axC, xticks=[0, 1500, 3000], xlim=(0, 3000), yticks=[0, 350, 700], ylim=(-45, 1000))
axC.yaxis.tick_right(); axC.yaxis.set_label_position("right")
axC.set_xlabel("practice step $t$", labelpad=1); axC.set_ylabel("$T$, ms", labelpad=2.5)
SUB = [S.text(BX + BW / 2, BY + BH + 0.05, "at one step, straight", fontsize=6.4, ha="center", va="bottom"),
       S.text(BX2 + BW2 / 2, BY + BH + 0.05, "over practice, curved", fontsize=6.4, ha="center", va="bottom")]

# ---------------------------------------------------------------- B, the two networks, one in each plot
# The even ring of step 0 sits above the line of step 0 at one step, and the coupling of the simulated performer of
# Figures 6 to 9 at step 1,000 sits above the 7-bit course at step 1,000. Each is drawn in the reading of Figure 2: node
# area by share, node colour by log2(p_i N), link width and shade by the number of links between the pair. The ring is a
# closed chain; in the performer's network its three hubs sit at the centre and the other elements are placed by a
# spring layout (networkx), its elements without links set as small open circles below it.
A1000 = sf.run(N=N, n_start=N, steps=1000, lam=0.5, t_elim=200, seed=int(SC["seed"]), return_adjacency=True)[3]
assert np.allclose(A1000.sum(1) / A1000.sum(), SC["p1000"])          # the shares of Section S4.4
H1000 = float(-(lambda q: (q * np.log2(q)).sum())(SC["p1000"][SC["p1000"] > 0]))
LINK, RIM = "#5E6E75", "#263238"
W_NET = 0.46


def layout(A):
    """Positions in [-1, 1]: a closed chain for the ring, otherwise the three hubs held apart at the centre, since they
    are joined to nearly every element, and the other elements placed around them by a spring layout."""
    import networkx as nx
    k = A.sum(1); on = [i for i in range(len(A)) if k[i] > 0]
    if np.all(k == k[0]) and np.all((A > 0).sum(1) == 2):
        G = nx.from_numpy_array(A); order = [0]
        while len(order) < len(A):
            order.append(next(v for v in G.neighbors(order[-1]) if v not in order))
        xy = loop_layout(len(A), 8)
        return {v: xy[i] / np.abs(xy).max() for i, v in enumerate(order)}
    G = nx.Graph(); G.add_nodes_from(on)
    for i, j in np.argwhere(np.triu(A, 1) > 0):
        G.add_edge(int(i), int(j))
    hubs = [int(i) for i in np.argsort(-k, kind="stable")[:3]]          # the three largest shares, at the centre
    init = {h: 0.45 * np.array([np.cos(a), np.sin(a)]) for h, a in zip(hubs, np.deg2rad([90, 210, 330]))}
    rng = np.random.default_rng(21)
    for i in on:
        if i not in init:
            init[i] = rng.uniform(-1, 1, 2)
    pos = nx.spring_layout(G, pos=init, fixed=hubs, k=0.6, iterations=400, seed=21)
    P = np.array([pos[i] for i in on]) - np.mean([pos[h] for h in hubs], 0)
    r = np.hypot(*P.T)[:, None]; P = P * (r / r.max()) ** -0.45         # draw the outer elements closer
    P /= np.abs(P).max()
    return {i: P[j] for j, i in enumerate(on)}


def draw_network(a_, A):
    n = len(A); k = A.sum(1); p = k / k.sum(); pos = layout(A)
    a_.set_xlim(-1.3, 1.3); a_.set_ylim(-1.4, 1.2); a_.set_aspect("equal"); a_.axis("off"); a_.patch.set_alpha(0.0)
    pairs = [(int(i), int(j)) for i, j in np.argwhere(np.triu(A, 1) > 0)]
    mm = max(A[i, j] for i, j in pairs)
    for i, j in sorted(pairs, key=lambda e: A[e]):
        f = np.sqrt(A[i, j] / mm)
        a_.plot(*zip(pos[i], pos[j]), color=LINK, lw=0.3 + 1.1 * f, alpha=0.40 + 0.55 * f, solid_capstyle="round",
                zorder=2)
    lone = [i for i in range(n) if k[i] == 0]
    for j, i in enumerate(lone):
        a_.scatter(0.30 * (j - (len(lone) - 1) / 2), -1.25, s=2.6, color="white", edgecolors=GREY, linewidths=0.35,
                   zorder=5)
    for i in sorted(pos, key=lambda z: p[z]):
        a_.scatter(*pos[i], s=2.2 + 95 * p[i], color=[share_colour(np.log2(p[i] * n))], edgecolors=RIM,
                   linewidths=0.4, zorder=5)


def network(axp, xd, yd, A, l1, l2, mk, target, w=W_NET):
    """A network of side w inches centred on (xd, yd) of the plot axp, its label of two lines just above it, the
    first led by the marker of its point in the plot, and a dotted arrow from the network to that point, as in Figures
    4B and 5B."""
    a_ = inset_at(S, axp, xd, yd, w); draw_network(a_, A)
    t1 = a_.text(0.5, 1.22, l1, transform=a_.transAxes, fontsize=6.2, ha="center", va="center")
    a_.text(0.5, 1.05, l2, transform=a_.transAxes, fontsize=6.2, ha="center", va="center")
    S.fig.canvas.draw()
    bb = t1.get_window_extent(S.fig.canvas.get_renderer()).transformed(a_.transAxes.inverted())
    a_.plot(bb.x0 - 0.09, 1.22, mk, ms=3.2, color="white" if mk == "o" else INK, mec=INK,
            mew=0.8 if mk == "o" else 0.3, transform=a_.transAxes, clip_on=False)
    axp.plot(*target, mk, ms=3.4 if mk == "o" else 3.2, color="white" if mk == "o" else INK,
             mec=INK if mk == "o" else "white", mew=0.8 if mk == "o" else 0.3, zorder=7)
    X0, Y0_ = axp.transData.transform((xd, yd)); X1, Y1_ = axp.transData.transform(target)
    r = 0.40 * w * S.fig.dpi; ang = np.arctan2(Y1_ - Y0_, X1 - X0)
    start = axp.transData.inverted().transform((X0 + r * np.cos(ang), Y0_ + r * np.sin(ang)))
    axp.annotate("", xy=target, xytext=start, zorder=6.5, annotation_clip=False,
                 arrowprops=dict(arrowstyle="-|>,head_length=0.32,head_width=0.16", color=EDGE, lw=0.6,
                                 ls=(0, (1.2, 1.4)), shrinkA=0, shrinkB=2.2))
    return a_


RING = np.zeros((N, N))
for i_ in range(N):
    RING[i_, (i_ + 1) % N] += 1; RING[(i_ + 1) % N, i_] += 1
T0 = (0.95, A_MS + BETA * H[0] * 0.95)                               # step 0 on its line
T1 = (1000, A_MS + BETA * H[1000] * 7)                               # step 1,000 on the 7-bit course
NETS = [network(axB, 1.55, 715, RING, "step 0", "$H$ = {:.2f} bits".format(CEIL), "o", T0),
        network(axC, 1650, 700, A1000, "step 1,000", "$H$ = {:.2f} bits".format(H1000), "D", T1, w=0.52)]


S.head("A", "A systems model of the dynamic law", AX0 + AWID / 2)
S.head("B", "One law through practice", (BX + BX2 + BW2) / 2)
out = ROOT / "figures" / "fig11_system_model.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [ax], "B": [axB, axC] + NETS + SUB}, str(ROOT / "figures" / "panels" / "fig11"))
print(f"H(0) {H[0]:.3f}, H(1000) {H[1000]:.3f}, H(3000) {H[3000]:.3f}; network of step 1,000: seed {int(SC['seed'])}, "
      f"H = {H1000:.4f} bits, {int((A1000.sum(1) > 0).sum())} elements with links; time saved by step 3,000: "
      f"{save[1]:.1f} ms at 1 bit, {save[7]:.1f} ms at 7 bits")
print("written", out)
