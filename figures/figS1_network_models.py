"""Figure S1 of the Supplemental Material: the network models of Section S2.1 and how far averaged information lies
below its ceiling.

A, the random, small-world and scale-free rules of Table S3 drawn at 60 elements with a mean degree of 4, built by
code/network_models.py from one seed, node area and grey by degree. The first two place every element at the start and
draw each link with a probability that ignores the links present, and the third grows by preferential attachment.
B, the share of elements with degree k or more at 1,000 elements, pooled over the 20 runs of Table S3 (the same
seeds), with k^-2, the cumulative form of P(k) ~ k^-3, for reference. A and B are drawn in grey.
For any degree sequence, log2 N - H = D, the relative entropy between the degree distribution seen from a link end,
q(k) = k P(k) / <k>, and the degree distribution P(k) itself (Equation S2). C, D against the exponent gamma of a
power law P(k) ~ k^-gamma for k >= m, in the continuum (Equation S3), with the exact gap of preferential growth,
P(k) = 2m(m + 1) / [k(k + 1)(k + 2)], for m = 1, 2 and 3, and the gap of a random network of mean degree 2.
D, the gap against the number of elements N for networks grown by preferential attachment with m = 1 (networkx), for
networks grown with fitness, 15 percent of elements at fitness 1.5 and the rest at 1, and for the four ensembles of
20 elements in Figure 2, made by the generators of Park (2026).
Usage: python -I figures/figS1_network_models.py <clone of github.com/pcw8531/control-entropy-coordinative-structures at 0660a649f5e0>
"""
import subprocess, sys
from pathlib import Path
import numpy as np
import matplotlib as mpl
import networkx as nx
from scipy.stats import poisson

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import WEIGHT, Sheet, INK, GREY, EDGE, RULE, TONE, CLEAR, tb_axes, save_panels  # noqa: E402
sys.path.insert(0, str(ROOT / "code"))
import network_models as nm  # noqa: E402  the rules of Table S3

DEP = Path(sys.argv[1]).resolve()
head_ = subprocess.run(["git", "-C", str(DEP), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
assert head_.startswith("0660a649f5e0"), f"deposit is at {head_}"
sys.path.insert(0, str(DEP / "code"))
from generators import make, SEEDS  # noqa: E402  Park (2026)

LN2 = np.log(2)


def gap(k):
    """D = log2 N - H for a degree sequence, elements with no link counted in N (Equation S2)."""
    k = np.asarray(k, float); N, K = len(k), k.sum(); kk = k[k > 0]
    return np.log2(N) - (np.log2(K) - (kk * np.log2(kk)).sum() / K)


def gap_power(g):
    """Continuum power law, k >= m, exponent g > 2 (Equation S3), in bits."""
    return (1 / (g - 2) - np.log((g - 1) / (g - 2))) / LN2


def gap_dist(k, P):
    P = P / P.sum(); mk = (k * P).sum(); q = k * P / mk
    return float((q[k > 0] * np.log2(k[k > 0] / mk)).sum())


kk = np.arange(1, 4_000_000, dtype=float)
BA = {m: gap_dist(kk[kk >= m], 2 * m * (m + 1) / (kk[kk >= m] * (kk[kk >= m] + 1) * (kk[kk >= m] + 2))) for m in (1, 2, 3)}
ER = gap_dist(np.arange(0, 80, dtype=float), poisson.pmf(np.arange(0, 80), 2.0))


def gap_fitness_continuum(etas, weights):
    """Continuum gap of growth with fitness (Equations S4 and S5): a mixture of power laws k^-(1 + C/eta), k >= 1."""
    from scipy.optimize import brentq
    etas, weights = np.asarray(etas, float), np.asarray(weights, float)
    C = brentq(lambda c: (weights * etas / (c - etas)).sum() - 1, etas.max() * (1 + 1e-9), 100.0)
    al = C / etas; Ek = (weights * al / (al - 1)).sum(); Ekl = (weights * al / (al - 1) ** 2).sum()
    return C, (Ekl / Ek - np.log(Ek)) / LN2


C_FIT, D_FIT = gap_fitness_continuum([1.0, 1.5], [0.85, 0.15])


def grow_fitness(N, share, eta_hi, rng):
    """Growth with m = 1 and two fitness classes. The class of the target is drawn by its total weight eta * K, and the
    element within the class by its degree, from a list holding one entry per link end."""
    hi = rng.random(N) < share; ends = {True: [], False: []}; K = {True: 0, False: 0}; k = np.zeros(N, int)
    for i in (0, 1):
        k[i] = 1; ends[bool(hi[i])].append(i); K[bool(hi[i])] += 1
    for t in range(2, N):
        w_hi = eta_hi * K[True]; c = rng.random() < w_hi / (w_hi + K[False])
        j = ends[c][rng.integers(len(ends[c]))]
        k[j] += 1; ends[bool(hi[j])].append(j); K[bool(hi[j])] += 1
        k[t] = 1; ends[bool(hi[t])].append(t); K[bool(hi[t])] += 1
    return k


rng = np.random.default_rng(20261106)
NS = [100, 1000, 10_000, 100_000]
sim_ba = {N: [gap([d for _, d in nx.barabasi_albert_graph(N, 1, seed=int(rng.integers(2**31))).degree()])
              for _ in range(20 if N < 100_000 else 5)] for N in NS}
sim_fit = {N: [gap(grow_fitness(N, 0.15, 1.5, rng)) for _ in range(20 if N < 100_000 else 5)] for N in NS}
ens = {kind: float(np.mean([gap([d for _, d in make(kind, 20, 2, s).degree()]) for s in SEEDS]))
       for kind in ("regular", "small-world", "random", "scale-free")}

FH = 4.35
S = Sheet(6.5, FH)
axC = S.ax(0.55, 0.42, 2.35, 1.50)
axD = S.ax(3.85, 0.42, 2.35, 1.50)

# ---------------------------------------------------------------- A, three rules drawn at 60 elements, in grey
SHOW = ("random", "small world", "preferential")
NAME = {"random": "random", "small world": "small world", "preferential": "scale-free"}
RULE_TXT = {"random": "static, each pair joined\nwith one probability",
            "small world": "static, each link of a ring\nrewired with probability 0.1",
            "preferential": "dynamic, growth with each link\ndrawn by $k_j/\\Sigma k$"}
nets60 = nm.build(nm.seeds()[0], n=60)
kmax60 = max(max(d for _, d in nets60[m][0].degree()) for m in SHOW)
axA = S.ax(0.05, 2.30, 4.55, 1.70); axA.set_xlim(0, 3); axA.set_ylim(-0.24, 1.02); axA.axis("off")
for c, name in enumerate(SHOW):
    G = nets60[name][0]
    pos = nx.circular_layout(G) if name == "small world" else nx.spring_layout(G, seed=3, k=0.22, iterations=500)
    xy = np.array([pos[v] for v in G.nodes()]); xy = (xy - xy.min(0)) / (xy.max(0) - xy.min(0)).max()
    xy = xy - (xy.max(0) + xy.min(0)) / 2
    P_ = {v: (c + 0.5 + 0.60 * xy[i, 0], 0.46 + 0.72 * xy[i, 1]) for i, v in enumerate(G.nodes())}
    for u, v in G.edges():
        axA.plot([P_[u][0], P_[v][0]], [P_[u][1], P_[v][1]], color="#A6A6A6", lw=0.45, zorder=1, solid_capstyle="round")
    for v, d in G.degree():
        g_ = 0.92 - 0.82 * d / kmax60
        axA.scatter(*P_[v], s=3 + 2.2 * d, color=(g_, g_, g_), edgecolor=INK, linewidth=0.3, zorder=2)
    axA.text(c + 0.5, 0.92, NAME[name], fontsize=7.4, ha="center", va="bottom", fontweight=WEIGHT)
    axA.text(c + 0.5, 0.02, RULE_TXT[name], fontsize=6.2, ha="center", va="top", linespacing=1.05)

# ---------------------------------------------------------------- B, degree distributions at 1,000 elements, in grey
axB = S.ax(5.05, 2.62, 1.30, 1.42)
deg = {m: [] for m in SHOW}
for s_ in nm.seeds():
    nets = nm.build(s_)
    for m in SHOW:
        deg[m] += [d for _, d in nets[m][0].degree()]
MK = {"random": dict(marker="o", color="white", mec=INK, mew=0.7), "small world": dict(marker="s", color="#8C8C8C", mec="#8C8C8C", mew=0.0),
      "preferential": dict(marker="^", color=INK, mec=INK, mew=0.0)}
for m in SHOW:                                                       # the share of elements with degree k or more
    d_ = np.sort(np.asarray(deg[m])); k_ = np.unique(d_[d_ > 0])
    ccdf = 1.0 - np.searchsorted(d_, k_, side="left") / d_.size
    axB.plot(k_, ccdf, ls="none", ms=2.6, label=NAME[m], zorder=3, **MK[m])
kk_ = np.array([2.0, 200.0]); axB.plot(kk_, 4.0 * kk_ ** -2.0, color=INK, lw=0.7, ls=(0, (3.0, 1.6)))
axB.text(30, 4.0 * 30.0 ** -2.0 * 1.6, "$k^{-2}$", fontsize=6.4, ha="left", va="bottom")
axB.set_xscale("log"); axB.set_yscale("log"); axB.set_xlim(1, 300); axB.set_ylim(2e-5, 1.5)
axB.set_xticks([1, 10, 100]); axB.set_yticks([1e-4, 1e-2, 1])
axB.xaxis.set_minor_locator(mpl.ticker.NullLocator()); axB.yaxis.set_minor_locator(mpl.ticker.NullLocator())
for sp in ("left", "right"):
    axB.spines[sp].set_visible(False)
axB.tick_params(axis="both", length=2.5, width=0.6, labelsize=6.6)
axB.set_xlabel("degree $k$", labelpad=1); axB.set_ylabel("share with degree $\\geq k$", labelpad=1.5)
axB.legend(loc="lower left", fontsize=6.0, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.2, handlelength=1.0, frameon=False)

ax = axC
g = np.linspace(2.2, 5.0, 300)
ax.plot(g, gap_power(g), color=INK, lw=1.1)
ax.axhline(ER, color=RULE, lw=0.7, ls=(0, (2.0, 1.6)))
ax.text(4.95, ER + 0.05, f"random, mean degree 2, {ER:.2f}", fontsize=6.2, ha="right", va="bottom", color=GREY)
ax.plot(3, gap_power(3), "o", ms=3.6, color=INK, zorder=4)
ax.text(2.94, 0.30, f"continuum, {gap_power(3):.3f}", fontsize=6.2, ha="right", va="top")
for m, dy in zip((1, 2, 3), (0.0, 0.0, 0.0)):
    ax.plot(3, BA[m], "o", ms=3.4, color="white", mec=TONE["scale-free"], mew=0.8, zorder=5)
ax.text(3.08, BA[1] + 0.02, f"growth, $m$ = 1, 2, 3:  {BA[1]:.3f}, {BA[2]:.3f}, {BA[3]:.3f}", fontsize=6.2,
        ha="left", va="bottom", color=TONE["scale-free"])
ax.text(2.42, 2.55, r"$D(\gamma)$", fontsize=7.2, ha="left")
tb_axes(ax, xticks=[2, 3.5, 5], xlim=(2.0, 5.0), yticks=[0, 1.5, 3], ylim=(0, 3.1))
ax.set_xlabel(r"exponent $\gamma$ of $P(k) \propto k^{-\gamma}$", labelpad=1)
ax.set_ylabel(r"gap $D = \log_2 N - H$, bits", labelpad=1.5)

ax = axD
for d, mk, kw, lab in ((sim_ba, "o", dict(color=TONE["scale-free"], mec=EDGE, mew=0.4), "preferential growth"),
                       (sim_fit, "^", dict(color="white", mec=INK, mew=0.8), "growth with fitness")):
    ax.errorbar(NS, [np.mean(d[N]) for N in NS], yerr=[np.std(d[N]) for N in NS], fmt=mk, ms=3.4, elinewidth=0.6,
                capsize=1.5, ecolor=GREY, label=lab, zorder=4, **kw)
ax.axhline(BA[1], color=INK, lw=0.8, ls=(0, (3.0, 1.6)))
ax.text(1.1e5, BA[1] - 0.03, f"exact, {BA[1]:.3f}", fontsize=6.2, ha="right", va="top")
for kind, mk in (("regular", "s"), ("small-world", "D"), ("random", "v"), ("scale-free", "o")):
    ax.plot(20, ens[kind], mk, ms=3.4, color=TONE[kind], mec=EDGE, mew=0.5, zorder=5)
for kind, txt, y_ in (("scale-free", "scale-free", ens["scale-free"] + 0.015), ("random", "random", ens["random"] - 0.015),
                      ("small-world", "small-world and regular", 0.02)):
    ax.text(26, y_, txt, fontsize=6.0, ha="left", va="center", color=GREY, **CLEAR)
ax.text(17, 0.47, "Figure 2", fontsize=6.2, ha="left", va="bottom")
ax.set_xscale("log"); ax.set_xlim(10, 2e5); ax.set_xticks([1e1, 1e3, 1e5])
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
tb_axes(ax, yticks=[0, 0.5, 1.0], ylim=(-0.05, 1.05))
ax.set_xlabel("number of elements $N$", labelpad=1); ax.set_ylabel(r"gap $D$, bits", labelpad=1.5)
ax.legend(loc="upper left", fontsize=6.2, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.2, handlelength=1.0)

S.head("A", "Static and dynamic rules", 0.05 + 4.55 / 2)
S.head_over(axB, "B", "Degrees")
S.head("C", "The gap set by the degree distribution", 0.55 + 2.35 / 2, y=2.04)
S.head("D", "The gap against size", 3.85 + 2.35 / 2, y=2.04)
out = ROOT / "figures" / "figS1_network_models.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axA], "B": [axB], "C": [axC], "D": [axD]}, str(ROOT / "figures" / "panels" / "figS1"),
            rows=[["A", "B"], ["C", "D"]])
print("A: largest degree at 60 elements:", {NAME[m]: max(d for _, d in nets60[m][0].degree()) for m in SHOW})
print("B: largest degree pooled at 1,000 elements:", {NAME[m]: max(deg[m]) for m in SHOW})
print("continuum gap, gamma 2.25, 2.5, 3, 4:", [round(float(gap_power(x)), 4) for x in (2.25, 2.5, 3, 4)])
print("exact growth gap, m = 1, 2, 3:", {m: round(v, 4) for m, v in BA.items()}, "| random, mean degree 2:", round(ER, 4))
print("preferential growth:", {N: (round(float(np.mean(v)), 4), round(float(np.std(v)), 4)) for N, v in sim_ba.items()})
print("growth with fitness:", {N: (round(float(np.mean(v)), 4), round(float(np.std(v)), 4)) for N, v in sim_fit.items()})
print(f"growth with fitness, continuum: C = {C_FIT:.4f}, tail exponent {1 + C_FIT / 1.5:.3f}, gap {D_FIT:.4f}")
print("Figure 2 ensembles, N = 20:", {k: round(v, 4) for k, v in ens.items()})
print("written", out)
