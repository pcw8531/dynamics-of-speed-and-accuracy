"""Where coupling goes during practice: the shares of three groups of elements, for the ternary panels of Figures 4B
and 5B of the manuscript (drawn by figures/fig4_minimal_model.py and figures/fig5_full_model.py).

The groups are read from the degrees at every step of every run:
  fit elements   the elements that practice favours, g_i > 0 (elements 0, 1 and 2, and in the late-recruit run also
                 the recruited element 19)
  chance hubs    the three strongest of the other elements at that step
  peripheral     all remaining elements, including any element not yet recruited (share zero)
Each group's share is its link ends divided by all link ends, so the three shares add to one. The panels draw the
fitness-weighted runs of Figure 4 and the fixed condition of Figure 5 (code/sim_flows.py takes them as its paths from
the even ring); the other runs are kept for the checks below.

The runs are those of the paper, with the same seeds and the same random draws, so the averaged information of every
run is unchanged. This script only records the degrees that the original scripts do not keep, and it checks that the
mean H(t) it obtains equals the stored mean H(t) of sim_out/sim_full_results.npz to machine precision.
  Figure 4 (code/sim_practice.py): random, preferential and fitness-weighted attachment, 300 runs to step 2,000, and
           the late recruit joining at step 50 with eta = 3, 200 runs to step 3,000.
  Figure 5 (code/sim_full.py): the fixed condition, the recruitment condition and the fixed set without elimination,
           300 runs to step 3,000.
Usage: python -I code/sim_composition.py sim_out   (about three minutes; writes sim_out/sim_composition.npz)
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import sim_practice as sp                                             # noqa: E402
import sim_full as sf                                                 # noqa: E402

FIT = [0, 1, 2]


def groups(K, fit):
    """Shares of (fit elements, chance hubs, peripheral elements) for a degree history K of shape (steps + 1, N)."""
    K = np.asarray(K, float)
    tot = K.sum(1)
    others = np.setdiff1d(np.arange(K.shape[1]), fit)
    f = K[:, fit].sum(1)
    top = np.sort(K[:, others], axis=1)[:, -3:].sum(1)
    return np.stack([f / tot, top / tot, (tot - f - top) / tot], 1)


def run_full_tracked(N=20, n_start=20, steps=3000, lam=0.5, g=None, t_elim=None, recruits=(), m=None, seed=0):
    """sim_full.run with the same random draws, keeping the degrees at every step."""
    rng = np.random.default_rng(seed)
    if g is None:
        g = np.zeros(N); g[:3] = 1.0
    eta = 1.0 + lam * np.asarray(g, float)
    A = np.zeros((N, N))
    for i in range(n_start):
        j = (i + 1) % n_start
        A[i, j] += 1; A[j, i] += 1
    k = A.sum(1)
    sched = {s: (idx, gg) for s, idx, gg in recruits}
    H = np.empty(steps + 1); H[0] = sf.H_bits(k)
    Kh = np.empty((steps + 1, N)); Kh[0] = k
    for t in range(1, steps + 1):
        if t in sched:
            idx, gg = sched[t]
            eta[idx] = 1.0 + lam * gg
            mm = m if m is not None else max(1, int(round(k.sum() / (k > 0).sum())))
            w = np.where(k > 0, eta * k, 0.0)
            for _ in range(mm):
                j = sf.draw(w, rng.random())
                A[idx, j] += 1; A[j, idx] += 1
            k = A.sum(1); H[t] = sf.H_bits(k); Kh[t] = k
            continue
        w = np.where(k > 0, eta * k, 0.0)
        i = sf.draw(w, rng.random()); w2 = w.copy(); w2[i] = 0.0; j = sf.draw(w2, rng.random())
        A[i, j] += 1; A[j, i] += 1; k = A.sum(1)
        if t_elim is not None and t >= t_elim:
            inv = np.where(k > 0, 1.0 / eta, 0.0)
            a = sf.draw(k * inv, rng.random())
            b = sf.draw(A[a], rng.random())
            A[a, b] -= 1; A[b, a] -= 1; k = A.sum(1)
        H[t] = sf.H_bits(k); Kh[t] = k
    return H, Kh


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    res, lines = {}, []

    # Figure 4, minimal model
    for rule in ("random", "pa", "fitness"):
        G, Hs = [], []
        for s in range(300):
            H, K = sp.run(rule, seed=s)
            G.append(groups(K, FIT)); Hs.append(H)
        res["m4_" + rule] = np.mean(G, 0); res["m4H_" + rule] = np.mean(Hs, 0)
        c = res["m4_" + rule]
        lines.append(f"Fig 4 {rule:8s} shares (fit, hubs, rest) at 0 {np.round(c[0], 3)}, 2000 {np.round(c[-1], 3)}; "
                     f"mean H(2000) {res['m4H_' + rule][-1]:.3f}")
    G = []
    for s in range(200):
        H, K = sp.run("fitness", steps=3000, seed=s, late=(19, 50, 4.0))
        G.append(groups(K, FIT + [19]))
    res["m4_late"] = np.mean(G, 0)
    c = res["m4_late"]
    lines.append(f"Fig 4 late     shares at 0 {np.round(c[0], 3)}, 50 {np.round(c[50], 3)}, 3000 {np.round(c[-1], 3)}")

    # Figure 5, full model, checked against the stored mean H(t)
    stored = np.load(ROOT / "sim_out" / "sim_full_results.npz")
    reg = sf.regimes()
    conds = {"fixed": reg["fixed"], "freeing": reg["freeing"],
             "fixed_no_elimination": dict(n_start=20, t_elim=None, recruits=())}
    for name, kw in conds.items():
        G, Hs = [], []
        for s in range(300):
            H, K = run_full_tracked(steps=3000, seed=s, **kw)
            G.append(groups(K, FIT)); Hs.append(H)
        Hm = np.mean(Hs, 0)
        dev = float(np.max(np.abs(Hm - stored[name])))
        res["f5_" + name] = np.mean(G, 0)
        c = res["f5_" + name]
        lines.append(f"Fig 5 {name:21s} shares at 0 {np.round(c[0], 3)}, 250 {np.round(c[250], 3)}, "
                     f"1000 {np.round(c[1000], 3)}, 3000 {np.round(c[-1], 3)}; "
                     f"largest gap to the stored mean H(t) {dev:.1e} bits")
        assert dev < 1e-9, "tracked runs differ from the stored runs"
    np.savez(out / "sim_composition.npz", **res)
    (out / "sim_composition_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
