"""The network models of Supplemental Section S2.1 and Table S3: which rules change a network by its own state.

Five rules build networks of N = 1,000 elements with a mean degree of 4, 20 runs each.
  ring lattice      every element linked to its 2 nearest neighbours on each side (Watts & Strogatz, 1998, p = 0)
  random            every pair linked with probability 4 / (N - 1) (Erdos & Renyi, 1959)
  small world       the ring lattice with each link rewired with probability 0.1 (Watts & Strogatz, 1998)
  preferential      growth from a triangle, each new element bringing m = 2 links whose other ends are drawn with
                    probability k_j / sum k (Barabasi & Albert, 1999; Equation 11)
  with fitness      the same growth with probability eta_j k_j / sum eta k, 15 percent of elements at fitness 1.5 and
                    the rest at 1, as in Figure S1 (Bianconi & Barabasi, 2001; Equation 12)
The first three place every element at the start and draw each link with a probability that does not depend on the
links already present. The last two draw each link by the degrees already present. For each network the script
records the largest degree divided by the mean degree, the mean clustering coefficient, the mean shortest path in the
largest connected component, and the gap of averaged information below its ceiling, log2 N - H, with H computed from
degree shares over all N elements (Equations 6, S1 and S2). For the two growing rules it also records the step at which
the element of largest degree joined.
Writes sim_out/network_models_result.txt and sim_out/network_models.csv.
The figure script of Figure S1 imports grow, build and seeds to draw the same rules.
Run from the repository root: python code/network_models.py sim_out
"""
import csv
import sys
from pathlib import Path
import numpy as np
import networkx as nx

N, K_MEAN, RUNS, SEED = 1000, 4, 20, 20261109
M = K_MEAN // 2


def gap(k):
    """log2 N - H for a degree sequence, every element counted in N (Equation S2)."""
    k = np.asarray(k, float); KS = k.sum(); kk = k[k > 0]
    return np.log2(len(k)) - (np.log2(KS) - (kk * np.log2(kk)).sum() / KS)


def grow(eta, rng, n=N, m=M):
    """Growth from a triangle, each new element bringing m links that attach with probability eta_j k_j / sum eta k
    to m different elements. Element t joins at step t, so its number is the step at which it joined."""
    G = nx.complete_graph(3); k = np.zeros(n); k[:3] = 2
    for t in range(3, n):
        w = eta[:t] * k[:t]
        targets = rng.choice(t, size=m, replace=False, p=w / w.sum())
        for j in targets:
            G.add_edge(t, int(j)); k[j] += 1
        k[t] = m
    return G


def build(seed, n=N):
    """The five networks of one run, built from one seed, as (graph, whether it grows)."""
    nets = {"ring lattice": (nx.watts_strogatz_graph(n, K_MEAN, 0.0, seed=seed), None),
            "random": (nx.gnp_random_graph(n, K_MEAN / (n - 1), seed=seed), None),
            "small world": (nx.watts_strogatz_graph(n, K_MEAN, 0.1, seed=seed), None),
            "preferential": (grow(np.ones(n), np.random.default_rng(seed), n), True)}
    eta = np.where(np.random.default_rng(seed + 1).random(n) < 0.15, 1.5, 1.0)
    nets["with fitness"] = (grow(eta, np.random.default_rng(seed), n), True)
    return nets


def seeds():
    """The seeds of the 20 runs of Table S3."""
    rng = np.random.default_rng(SEED)
    return [int(rng.integers(2**31)) for _ in range(RUNS)]


def stats(G, joined=None):
    deg = np.array([d for _, d in sorted(G.degree())], float)
    big = G.subgraph(max(nx.connected_components(G), key=len))
    row = {"k_mean": deg.mean(), "kmax_over_kmean": deg.max() / deg.mean(),
           "clustering": nx.average_clustering(G), "path_length": nx.average_shortest_path_length(big),
           "gap_bits": gap(deg)}
    if joined is not None:
        row["hub_joined_step"] = int(np.argmax(deg))                 # elements are numbered by the step they joined
    return row


MODELS = ["ring lattice", "random", "small world", "preferential", "with fitness"]


def main(out):
    rows = []
    for run, s in enumerate(seeds()):
        nets = build(s)
        for name in MODELS:
            G, j = nets[name]
            assert G.number_of_nodes() == N
            rows.append({"model": name, "run": run, "seed": s, **stats(G, j)})

    out.mkdir(exist_ok=True)
    keys = ["model", "run", "seed", "k_mean", "kmax_over_kmean", "clustering", "path_length", "gap_bits", "hub_joined_step"]
    with open(out / "network_models.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})

    lines = [f"Network models, N = {N}, mean degree {K_MEAN}, {RUNS} runs each, seed {SEED}. Means over runs, "
             "with the smallest and largest run in brackets.", ""]
    for name in MODELS:
        rr = [r for r in rows if r["model"] == name]

        def fmt(key, d):
            v = np.array([r[key] for r in rr], float)
            return f"{v.mean():.{d}f} [{v.min():.{d}f}, {v.max():.{d}f}]"
        line = (f"{name:13s} mean degree {fmt('k_mean', 2)}, largest/mean {fmt('kmax_over_kmean', 1)}, "
                f"clustering {fmt('clustering', 3)}, path {fmt('path_length', 2)}, log2 N - H {fmt('gap_bits', 3)} bits")
        if "hub_joined_step" in rr[0]:
            line += (f", hub joined at step median {int(np.median([r['hub_joined_step'] for r in rr]))}, "
                     f"latest {max(r['hub_joined_step'] for r in rr)}")
        lines.append(line)
    (out / "network_models_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("sim_out"))
