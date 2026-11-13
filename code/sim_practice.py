"""Minimal practice model for the dynamic law (Figure 4 of the manuscript, drawn by
figures/fig4_minimal_model.py).

N elements start with even coupling (a ring, every element has 2 links), so H(0) = log2 N.
This is a fixed set of elements whose coupling grows; elements are not added (except the one late recruit).
Each practice step adds one link between two different elements, chosen by one of three rules:
  random  : every element equally likely (energy spent without direction)
  pa      : probability proportional to degree k_i (Equation 11)
  fitness : probability proportional to eta_i * k_i, eta_i = 1 + lam * g_i (Equations 12 and 13)
Optionally one element is recruited late, at step t_join, with fitness eta_late.
"""
import numpy as np


def H_bits(k):
    p = k[k > 0] / k[k > 0].sum()
    return float(-(p * np.log2(p)).sum())


def run(rule, N=20, steps=2000, lam=0.5, g=None, seed=0, late=None):
    """Return H(t) for t = 0..steps and the degree history (steps+1, N).
    g: array of length N with g_i >= 0 (default: g_i = 1 for elements 0, 1, 2, else 0).
    late: None or (index, t_join, g_late). The other N - 1 elements start as a ring among themselves
          (degree 2 each); at step t_join the recruit joins with one link to an element chosen by the rule."""
    rng = np.random.default_rng(seed)
    if g is None:
        g = np.zeros(N); g[:3] = 1.0
    eta = 1.0 + lam * np.asarray(g, float) if rule == "fitness" else np.ones(N)
    k = np.full(N, 2.0)
    active = np.ones(N, bool)
    if late is not None:
        idx, t_join, g_late = late
        active[idx] = False; k[idx] = 0.0
    H = np.empty(steps + 1); K = np.empty((steps + 1, N))
    H[0] = H_bits(k); K[0] = k
    u = rng.random((steps, 2))
    for t in range(1, steps + 1):
        if late is not None and t == late[1]:        # the recruit joins with one link
            w = active.astype(float) if rule == "random" else np.where(active, eta * k, 0.0)
            c = np.cumsum(w); j = min(np.searchsorted(c, u[t - 1, 0] * c[-1], side="right"), N - 1)
            active[late[0]] = True; k[late[0]] = 1.0; k[j] += 1
            eta[late[0]] = 1.0 + lam * late[2]
            H[t] = H_bits(k); K[t] = k
            continue
        w = active.astype(float) if rule == "random" else np.where(active, eta * k, 0.0)
        c = np.cumsum(w); i = min(np.searchsorted(c, u[t - 1, 0] * c[-1], side="right"), N - 1)
        w2 = w.copy(); w2[i] = 0.0
        c2 = np.cumsum(w2); j = min(np.searchsorted(c2, u[t - 1, 1] * c2[-1], side="right"), N - 1)
        k[i] += 1; k[j] += 1
        H[t] = H_bits(k); K[t] = k
    return H, K


if __name__ == "__main__":
    out = {}
    for rule in ["random", "pa", "fitness"]:
        Hs = np.array([run(rule, seed=s)[0] for s in range(300)])
        out["H_" + rule] = Hs.mean(0)
        print(f"{rule:8s} H(0)={Hs.mean(0)[0]:.3f}  H(100)={Hs.mean(0)[100]:.3f}  "
              f"H(500)={Hs.mean(0)[500]:.3f}  H(2000)={Hs.mean(0)[-1]:.3f}")
    # late recruit: element 19 joins at step 50 with g = 4 (eta = 3); 3000 steps, 200 runs
    for t_join in (50, 200):
        Ks = np.array([run("fitness", steps=3000, seed=s, late=(19, t_join, 4.0))[1] for s in range(200)])
        top = np.mean(Ks[:, -1, 19] > Ks[:, -1, :19].max(1))
        print(f"late element joining at {t_join}: strongest at step 3000 in {100*top:.0f}% of runs")
        if t_join == 50:
            out["K_late"] = Ks[:, :, 19].mean(0)
            out["K_target"] = Ks[:, :, :3].mean((0, 2))
            out["K_other"] = Ks[:, :, 3:19].mean((0, 2))
            out["top50"] = top
        else:
            out["top200"] = top
    # a late recruit catches up if practice goes on long enough
    Ks = np.array([run("fitness", steps=15000, seed=s, late=(19, 200, 4.0))[1] for s in range(100)])
    for T_end in (3000, 6000, 15000):
        top = np.mean(Ks[:, T_end, 19] > Ks[:, T_end, :19].max(1))
        print(f"late element joining at 200: strongest at step {T_end} in {100*top:.0f}% of 100 runs")
    np.savez("sim_results.npz", **out)
    a, b_beg, ID = 12.8, 94.7, 4.0
    beta = b_beg / np.log2(20)
    for r in ["random", "pa", "fitness"]:
        H = out["H_" + r]
        print(r, "MT(0)=%.1f MT(2000)=%.1f ms" % (a + beta * H[0] * ID, a + beta * H[-1] * ID))
    print("beta =", beta)
