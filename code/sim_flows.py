"""Where the coupling goes from other starts in the triangle of shares, for the middle panels of Figures 4B and 5B of
the manuscript (drawn by figures/fig4_minimal_model.py and figures/fig5_full_model.py).

The groups are those of code/sim_composition.py: the fit elements (0, 1 and 2), the three strongest of the other
elements at each step (chance hubs), and the peripheral elements, the 14 that remain. Each start gives each fit
element df link ends, each of three other elements dh and each of the remaining 14 elements dr, between 40 and 58 link
ends in all, close to the 40 of the even ring. In every start the peripheral elements hold the largest share of the
coupling at step 0, 0.48 to 0.77 (0.70 in the ring), so the starts lie toward the lower left corner of the triangle,
along the base toward the fit elements and up toward the chance hubs. The start (2, 2, 2) is the even ring of the
paper, whose paths are those of sim_out/sim_composition.npz. Every other start runs 300 times with seeds 0 to 299, the
same draws as the runs of the paper.

Figure 4 (minimal model, code/sim_practice.py): fitness-weighted attachment, lambda = 0.5, 2,000 steps. A start is a
           degree vector.
Figure 5 (full model, code/sim_full.py): the fixed set with selective elimination from step 200, 3,000 steps. A start
           is a multigraph with the given degrees, made by joining link ends one at a time without self-links
           (seed 7000 plus the index of the start).
The script first checks that its two runners reproduce the runs of the paper exactly from the even ring.
Usage: python -I code/sim_flows.py sim_out   (about two and a half minutes; writes sim_out/sim_flows.npz)
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import sim_practice as sp                                             # noqa: E402
import sim_full as sf                                                 # noqa: E402
from sim_composition import groups, FIT                               # noqa: E402

STARTS = [(2, 2, 2), (1, 2, 3), (1, 4, 2), (1, 6, 2), (3, 2, 2), (5, 2, 2), (8, 2, 2)]
RUNS = 300


def degrees(start, n_fit=3, n_hub=3, n_rest=14):
    df, dh, dr = start
    return np.array([df] * n_fit + [dh] * n_hub + [dr] * n_rest, float)


def run_min_from(k0, rule, steps=2000, lam=0.5, seed=0):
    """sim_practice.run with the starting degrees k0 instead of the ring; the same draws for the ring start."""
    rng = np.random.default_rng(seed)
    N = len(k0); g = np.zeros(N); g[:3] = 1.0
    eta = 1.0 + lam * g if rule == "fitness" else np.ones(N)
    k = np.array(k0, float)
    H = np.empty(steps + 1); K = np.empty((steps + 1, N)); H[0] = sp.H_bits(k); K[0] = k
    u = rng.random((steps, 2))
    for t in range(1, steps + 1):
        w = np.ones(N) if rule == "random" else eta * k
        c = np.cumsum(w); i = min(np.searchsorted(c, u[t - 1, 0] * c[-1], side="right"), N - 1)
        w2 = w.copy(); w2[i] = 0.0
        c2 = np.cumsum(w2); j = min(np.searchsorted(c2, u[t - 1, 1] * c2[-1], side="right"), N - 1)
        k[i] += 1; k[j] += 1
        H[t] = sp.H_bits(k); K[t] = k
    return H, K


def multigraph(deg, seed):
    """A multigraph without self-links whose degrees are deg (the total is made even on the last element). Link ends
    are paired one at a time: the element with most open ends is joined to another element drawn in proportion to its
    open ends, which never needs a self-link while no element holds more than half of all ends."""
    deg = np.array(deg, int)
    if deg.sum() % 2:
        deg[-1] += 1
    rng = np.random.default_rng(seed)
    left = deg.copy(); A = np.zeros((len(deg), len(deg)))
    while left.sum() > 0:
        i = int(np.argmax(left)); w = left.astype(float); w[i] = 0.0
        j = int(rng.choice(len(deg), p=w / w.sum()))
        A[i, j] += 1; A[j, i] += 1; left[i] -= 1; left[j] -= 1
    return A


def multigraph_ring(n):
    A = np.zeros((n, n))
    for i in range(n):
        j = (i + 1) % n; A[i, j] += 1; A[j, i] += 1
    return A


def run_full_from(A0, steps=3000, lam=0.5, t_elim=None, seed=0):
    """sim_full.run on the fixed set, starting from the multigraph A0 instead of the ring."""
    rng = np.random.default_rng(seed)
    N = A0.shape[0]; g = np.zeros(N); g[:3] = 1.0
    eta = 1.0 + lam * g
    A = A0.copy(); k = A.sum(1)
    H = np.empty(steps + 1); H[0] = sf.H_bits(k); Kh = np.empty((steps + 1, N)); Kh[0] = k
    for t in range(1, steps + 1):
        w = np.where(k > 0, eta * k, 0.0)
        i = sf.draw(w, rng.random()); w2 = w.copy(); w2[i] = 0.0; j = sf.draw(w2, rng.random())
        A[i, j] += 1; A[j, i] += 1; k = A.sum(1)
        if t_elim is not None and t >= t_elim:
            inv = np.where(k > 0, 1.0 / eta, 0.0)
            a = sf.draw(k * inv, rng.random()); b = sf.draw(A[a], rng.random())
            A[a, b] -= 1; A[b, a] -= 1; k = A.sum(1)
        H[t] = sf.H_bits(k); Kh[t] = k
    return H, Kh


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    comp = np.load(ROOT / "sim_out" / "sim_composition.npz")
    res, lines = {"starts": np.array(STARTS)}, []

    # from the even ring both runners give the runs of the paper
    for s in range(3):
        H0, K0 = sp.run("fitness", seed=s); H1, K1 = run_min_from(degrees((2, 2, 2)), "fitness", seed=s)
        assert np.array_equal(H0, H1) and np.array_equal(K0, K1)
    H0, _, _ = sf.run(steps=300, seed=1, n_start=20, t_elim=200)
    H1, _ = run_full_from(multigraph_ring(20), steps=300, t_elim=200, seed=1)
    assert np.array_equal(H0, H1)

    for si, st in enumerate(STARTS):
        if st == (2, 2, 2):
            res[f"m4_{si}"], res[f"m4H_{si}"] = comp["m4_fitness"], comp["m4H_fitness"]
        else:
            G, Hs = [], []
            for s in range(RUNS):
                H, K = run_min_from(degrees(st), "fitness", seed=s)
                G.append(groups(K, FIT)); Hs.append(H)
            res[f"m4_{si}"], res[f"m4H_{si}"] = np.mean(G, 0), np.mean(Hs, 0)
        c = res[f"m4_{si}"]
        lines.append(f"Fig 4 fitness-weighted, start {st}: shares {np.round(c[0], 3)} -> {np.round(c[1000], 3)} (1000) -> "
                     f"{np.round(c[-1], 3)} (2000), H {res[f'm4H_{si}'][0]:.3f} -> {res[f'm4H_{si}'][-1]:.3f}")
    for si, st in enumerate(STARTS):
        if st == (2, 2, 2):
            res[f"f5_{si}"] = comp["f5_fixed"]
            res[f"f5H_{si}"] = np.load(ROOT / "sim_out" / "sim_full_results.npz")["fixed"]
        else:
            A0 = multigraph(degrees(st), 7000 + si)
            G, Hs = [], []
            for s in range(RUNS):
                H, K = run_full_from(A0, seed=s, t_elim=200)
                G.append(groups(K, FIT)); Hs.append(H)
            res[f"f5_{si}"], res[f"f5H_{si}"] = np.mean(G, 0), np.mean(Hs, 0)
        c = res[f"f5_{si}"]
        lines.append(f"Fig 5 selective elimination, start {st}: shares {np.round(c[0], 3)} -> {np.round(c[1000], 3)} "
                     f"(1000) -> {np.round(c[-1], 3)} (3000), H {res[f'f5H_{si}'][0]:.3f} -> {res[f'f5H_{si}'][-1]:.3f}")
    np.savez(out / "sim_flows.npz", **res)
    (out / "sim_flows_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
