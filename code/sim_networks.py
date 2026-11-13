"""The networks above the triangles of Figures 4B and 5B: one run of each model from the even ring, drawn at the first
and the last step (figures/fig4_minimal_model.py, figures/fig5_full_model.py).

Each step of either model joins two elements by a link, so the links of a run form a multigraph that starts as the ring
of 20 elements. This script repeats the runs of the paper with the same seeds and draws, keeps the links, and checks
that the degrees (minimal model) and averaged information (full model) equal those of code/sim_practice.py and
code/sim_full.py. The run drawn is the one, among seeds 0 to 299, whose shares of fit elements, chance hubs and
peripheral elements at the last step lie closest to the mean of all 300 runs.

Minimal model: fitness-weighted attachment, lambda = 0.5, 2,000 steps (Figure 4B).
Full model: the fixed set with selective elimination from step 200, lambda = 0.5, 3,000 steps (Figure 5B).
Usage: python -I code/sim_networks.py sim_out   (about half a minute; writes sim_out/sim_networks.npz)
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import sim_practice as sp                                             # noqa: E402
import sim_full as sf                                                 # noqa: E402
from sim_composition import groups, FIT                               # noqa: E402

N = 20


def ring(n=N):
    A = np.zeros((n, n))
    for i in range(n):
        j = (i + 1) % n; A[i, j] += 1; A[j, i] += 1
    return A


def run_min_links(steps=2000, lam=0.5, seed=0):
    """sim_practice.run('fitness') with the same draws, keeping the links: the adjacency at the last step."""
    rng = np.random.default_rng(seed)
    g = np.zeros(N); g[:3] = 1.0; eta = 1.0 + lam * g
    A = ring(); k = A.sum(1)
    u = rng.random((steps, 2))
    for t in range(1, steps + 1):
        w = eta * k
        c = np.cumsum(w); i = min(np.searchsorted(c, u[t - 1, 0] * c[-1], side="right"), N - 1)
        w2 = w.copy(); w2[i] = 0.0
        c2 = np.cumsum(w2); j = min(np.searchsorted(c2, u[t - 1, 1] * c2[-1], side="right"), N - 1)
        A[i, j] += 1; A[j, i] += 1; k = A.sum(1)
    return A


def run_full_links(steps=3000, lam=0.5, t_elim=200, seed=0):
    """sim_full.run on the fixed set with the same draws, keeping the links: H(t) and the adjacency at the last step."""
    rng = np.random.default_rng(seed)
    g = np.zeros(N); g[:3] = 1.0; eta = 1.0 + lam * g
    A = ring(); k = A.sum(1)
    H = np.empty(steps + 1); H[0] = sf.H_bits(k)
    for t in range(1, steps + 1):
        w = np.where(k > 0, eta * k, 0.0)
        i = sf.draw(w, rng.random()); w2 = w.copy(); w2[i] = 0.0; j = sf.draw(w2, rng.random())
        A[i, j] += 1; A[j, i] += 1; k = A.sum(1)
        if t >= t_elim:
            inv = np.where(k > 0, 1.0 / eta, 0.0)
            a = sf.draw(k * inv, rng.random()); b = sf.draw(A[a], rng.random())
            A[a, b] -= 1; A[b, a] -= 1; k = A.sum(1)
        H[t] = sf.H_bits(k)
    return H, A


def shares(k):
    return groups(np.asarray(k, float)[None, :], FIT)[0]


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    res, lines = {}, []
    # minimal model: the degrees of every run, and the run closest to the mean composition at step 2,000
    K = np.array([sp.run("fitness", seed=s)[1][-1] for s in range(300)])
    C = np.array([shares(k) for k in K]); best = int(np.argmin(np.abs(C - C.mean(0)).sum(1)))
    A = run_min_links(seed=best)
    assert np.array_equal(A.sum(1), K[best])
    res.update(m4_seed=best, m4_A0=ring(), m4_A=A, m4_H=sp.H_bits(A.sum(1)))
    lines.append(f"Figure 4: seed {best}, shares at step 2,000 {np.round(C[best], 3)} (mean of 300 runs "
                 f"{np.round(C.mean(0), 3)}), H {sp.H_bits(A.sum(1)):.3f} bits, {int(A.sum() / 2)} links")
    # full model: the same for the fixed set with selective elimination at step 3,000
    runs = [run_full_links(seed=s) for s in range(300)]
    for s in range(3):
        assert np.array_equal(runs[s][0], sf.run(N=20, n_start=20, steps=3000, lam=0.5, t_elim=200, seed=s)[0])
    C = np.array([shares(A_.sum(1)) for _, A_ in runs]); best = int(np.argmin(np.abs(C - C.mean(0)).sum(1)))
    H, A = runs[best]
    res.update(f5_seed=best, f5_A0=ring(), f5_A=A, f5_H=H[-1])
    lines.append(f"Figure 5: seed {best}, shares at step 3,000 {np.round(C[best], 3)} (mean of 300 runs "
                 f"{np.round(C.mean(0), 3)}), H {H[-1]:.3f} bits, {int(A.sum() / 2)} links, "
                 f"{int((A.sum(1) == 0).sum())} elements without links")
    np.savez(out / "sim_networks.npz", **res)
    (out / "sim_networks_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
