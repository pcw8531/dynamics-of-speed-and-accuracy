"""Full practice model of the dynamic law: formation, selective elimination and recruitment.

The analytic result checked here (Supplemental Material, Section S3.2):
adding one unit of coupling to element i changes H, to first order, by (I_i - H) / K,
where I_i = -log2 p_i is the information of element i and K is the number of link ends.
Removing one unit changes H by -(I_i - H) / K. An element recruited with m links changes H by about
q [log2(1/q) + log2 e - H] + (m/K)(I_pi - H), with q = m / (K + 2m) its share at entry
and I_pi the mean of I_i over the elements its partners are drawn from.

Elements: N_max potential elements. An element with no links is inactive and holds share 0, so H is
always computed over the same fixed set of N_max elements (zero shares contribute nothing).

Events at each practice step
  formation   one link; first end drawn with weight eta_i k_i, second end by the same weights from the others
  elimination from step t_elim on, one link per step, so the number of links stays constant after the delay
              (Xu et al., 2009). Each link end carries removal weight 1/eta of its element; an end is drawn
              by that weight and its link is removed (selective stabilization: fit elements keep their links)
  recruitment at scheduled steps an inactive element joins with m links, partners drawn with weight eta_j k_j;
              by default m is the mean degree of the active elements, so the recruited element joins as a typical one

Usage: python -I code/sim_full.py <output folder>
"""
import sys
from pathlib import Path
import numpy as np

LOG2E = 1.0 / np.log(2.0)


def H_bits(k):
    k = k[k > 0]
    p = k / k.sum()
    return float(-(p * np.log2(p)).sum())


def info(k):
    """I_i = -log2 p_i for active elements, nan for inactive ones."""
    I = np.full(len(k), np.nan)
    a = k > 0
    I[a] = -np.log2(k[a] / k.sum())
    return I


# ---------- closed-form first-order changes (Section S3.2) ----------
def ends_formation(k, eta):
    """Expected number of ends added to each element by one formed link (sums to 2), exact for the rule."""
    w = np.where(k > 0, eta * k, 0.0)
    S = w.sum()
    P = np.outer(w / S, w) / (S - w)[:, None]       # P[i, j] = P(first end i) P(second end j | i)
    np.fill_diagonal(P, 0.0)
    return P.sum(1) + P.sum(0), P


def ends_elimination(A, k, eta):
    """Expected number of ends removed from each element by one eliminated link (sums to 2), exact for the rule.
    Link (i, j) is removed with probability proportional to its multiplicity times (1/eta_i + 1/eta_j)."""
    inv = np.where(k > 0, 1.0 / eta, 0.0)
    W = np.triu(A, 1) * (inv[:, None] + inv[None, :])
    P = W / W.sum()
    return P.sum(1) + P.sum(0), P


def first_order(k, e, sign=+1):
    """First-order change of H when e_i ends are added (sign +1) or removed (sign -1)."""
    K = k.sum(); H = H_bits(k); I = info(k)
    a = k > 0
    return sign * float(np.nansum(e[a] * (I[a] - H))) / K


def second_order(k, e, sign=+1):
    """Second-order term for ends that each change a degree by one (no self-links), n = 2 sign ends in all.
    In nats: -(1/2K) sum e_i/k_i + n^2/(2K^2) - (n/K^2) sum e_i (I_i - H). It is the chance spreading of the Polya urn."""
    K = k.sum(); H = H_bits(k); I = info(k); a = k > 0; n = 2.0 * sign
    quad = (-(np.sum(e[a] / k[a])) / (2 * K) + n ** 2 / (2 * K ** 2)) * LOG2E
    return float(quad - (n / K ** 2) * sign * np.nansum(e[a] * (I[a] - H)))


def recruit_first_order(k, eta, m):
    K = k.sum(); H = H_bits(k); I = info(k)
    w = np.where(k > 0, eta * k, 0.0); pi = w / w.sum()
    q = m / (K + 2 * m)
    return q * (np.log2(1 / q) + LOG2E - H) + (m / K) * (np.nansum(pi * I) - H)


# ---------- exact expectations by enumeration, for the check ----------
def exact_formation(k, eta):
    _, P = ends_formation(k, eta)
    H0 = H_bits(k); out = 0.0
    for i, j in zip(*np.nonzero(P)):
        kk = k.copy(); kk[i] += 1; kk[j] += 1
        out += P[i, j] * (H_bits(kk) - H0)
    return out


def exact_elimination(A, k, eta):
    _, P = ends_elimination(A, k, eta)
    H0 = H_bits(k); out = 0.0
    for i, j in zip(*np.nonzero(P)):
        kk = k.copy(); kk[i] -= 1; kk[j] -= 1
        out += P[i, j] * (H_bits(kk) - H0)
    return out


def exact_recruit(k, eta, m, new):
    """m = 1 or 2, partners drawn with replacement by weight eta_j k_j."""
    w = np.where(k > 0, eta * k, 0.0); pi = w / w.sum()
    H0 = H_bits(k); out = 0.0
    idx = np.nonzero(pi)[0]
    if m == 1:
        for j in idx:
            kk = k.copy(); kk[new] = 1; kk[j] += 1
            out += pi[j] * (H_bits(kk) - H0)
    else:
        for j in idx:
            for l in idx:
                kk = k.copy(); kk[new] = 2; kk[j] += 1; kk[l] += 1
                out += pi[j] * pi[l] * (H_bits(kk) - H0)
    return out


# ---------- the model ----------
def draw(w, u):
    c = np.cumsum(w)
    return min(int(np.searchsorted(c, u * c[-1], side="right")), len(w) - 1)


def run(N=20, n_start=20, steps=3000, lam=0.5, g=None, t_elim=None, recruits=(), m=None, seed=0,
        track_prediction=False, return_adjacency=False):
    """recruits: list of (step, element index, g of the recruited element). Elements n_start..N-1 start inactive.
    Returns H(t), degrees at the end, and (if asked) the first-order predicted change at every step. With
    return_adjacency=True the link counts between elements at the end are returned as a fourth value; the run itself
    is unchanged."""
    rng = np.random.default_rng(seed)
    if g is None:
        g = np.zeros(N); g[:3] = 1.0
    eta = 1.0 + lam * np.asarray(g, float)
    A = np.zeros((N, N))
    for i in range(n_start):                        # ring among the starting elements
        j = (i + 1) % n_start
        A[i, j] += 1; A[j, i] += 1
    k = A.sum(1)
    sched = {s: (idx, gg) for s, idx, gg in recruits}
    H = np.empty(steps + 1); H[0] = H_bits(k)
    pred = np.zeros(steps + 1) if track_prediction else None
    for t in range(1, steps + 1):
        if t in sched:                              # recruitment replaces formation at this step
            idx, gg = sched[t]
            eta[idx] = 1.0 + lam * gg
            mm = m if m is not None else max(1, int(round(k.sum() / (k > 0).sum())))
            if track_prediction:
                pred[t] = recruit_first_order(k, eta, mm)
            w = np.where(k > 0, eta * k, 0.0)
            for _ in range(mm):
                j = draw(w, rng.random())
                A[idx, j] += 1; A[j, idx] += 1
            k = A.sum(1); H[t] = H_bits(k)
            continue
        if track_prediction:
            e, _ = ends_formation(k, eta)
            pred[t] = first_order(k, e, +1) + second_order(k, e, +1)
        w = np.where(k > 0, eta * k, 0.0)
        i = draw(w, rng.random()); w2 = w.copy(); w2[i] = 0.0; j = draw(w2, rng.random())
        A[i, j] += 1; A[j, i] += 1; k = A.sum(1)
        if t_elim is not None and t >= t_elim:
            if track_prediction:
                e, _ = ends_elimination(A, k, eta)
                pred[t] += first_order(k, e, -1) + second_order(k, e, -1)
            inv = np.where(k > 0, 1.0 / eta, 0.0)
            a = draw(k * inv, rng.random())         # an end, by weight 1/eta of its element
            b = draw(A[a], rng.random())            # its link, by multiplicity
            A[a, b] -= 1; A[b, a] -= 1; k = A.sum(1)
        H[t] = H_bits(k)
    if return_adjacency:
        return H, k, pred, A
    return H, k, pred


def regimes(N=20):
    """Practice on a body of N_max = 20 elements, fitness lam = 0.5 with g = 1 on elements 0 to 2.
    fixed      all 20 elements active from the start, formation, and selective elimination from step 200
    freeing    10 elements active at first (a frozen beginner), the other 10 recruited one by one at steps
               25, 50, ..., 250 with g = 0, each joining with the mean degree, then the same formation and elimination"""
    return {
        "fixed": dict(n_start=N, t_elim=200, recruits=()),
        "freeing": dict(n_start=10, t_elim=200, recruits=[(25 * s, 9 + s, 0.0) for s in range(1, 11)]),
    }


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    lines = []

    # 1. single-event check at random states: closed form against exact enumeration
    rng = np.random.default_rng(20261105)
    errs = {"formation": [], "elimination": [], "recruitment": []}
    for trial in range(200):
        N = int(rng.integers(8, 41))
        H_, kk, _ = run(N=N, n_start=N - 1, steps=int(rng.integers(20, 400)), lam=float(rng.uniform(0, 2)),
                        g=rng.integers(0, 3, N).astype(float), t_elim=None, seed=trial)
        eta = 1.0 + float(rng.uniform(0, 2)) * rng.integers(0, 3, N)
        # rebuild an adjacency consistent with kk for the elimination check: use a fresh run that keeps A
        e, _ = ends_formation(kk, eta)
        errs["formation"].append((first_order(kk, e, +1) + second_order(kk, e, +1), exact_formation(kk, eta),
                                  first_order(kk, e, +1)))
        m = int(rng.integers(1, 3)); new = N - 1
        f = recruit_first_order(kk, eta, m)
        errs["recruitment"].append((f, exact_recruit(kk, eta, m, new), f))
    # elimination needs the adjacency, so run the model directly here
    for trial in range(200):
        N = int(rng.integers(8, 41)); eta = 1.0 + float(rng.uniform(0, 2)) * rng.integers(0, 3, N)
        A = np.zeros((N, N))
        for i in range(N):
            j = (i + 1) % N; A[i, j] += 1; A[j, i] += 1
        for _ in range(int(rng.integers(20, 400))):
            k = A.sum(1); w = eta * k
            i = draw(w, rng.random()); w2 = w.copy(); w2[i] = 0; j = draw(w2, rng.random())
            A[i, j] += 1; A[j, i] += 1
        k = A.sum(1); e, _ = ends_elimination(A, k, eta)
        errs["elimination"].append((first_order(k, e, -1) + second_order(k, e, -1), exact_elimination(A, k, eta),
                                first_order(k, e, -1)))
    lines.append("1. Closed form against exact enumeration, 200 random states each "
                 "(formation and elimination: first plus second order; recruitment: first order)")
    for name, v in errs.items():
        v = np.array(v); r = np.corrcoef(v[:, 0], v[:, 1])[0, 1]
        slope = np.polyfit(v[:, 1], v[:, 0], 1)[0]
        gap = np.median(np.abs(v[:, 0] - v[:, 1])) / np.median(np.abs(v[:, 1]))
        sign = np.mean(np.sign(v[:, 0]) == np.sign(v[:, 1]))
        r1 = np.corrcoef(v[:, 2], v[:, 1])[0, 1]
        lines.append(f"   {name:9s} r = {r:.4f}, slope {slope:.3f}, median gap {100 * gap:.1f}% of the median change, "
                     f"same sign {100 * sign:.0f}%; first order alone r = {r1:.4f}")

    # 2. three regimes, 300 runs each
    reg = regimes()
    curves = {}
    for name, kw in reg.items():
        Hs = np.array([run(steps=3000, seed=s, **kw)[0] for s in range(300)])
        curves[name] = Hs.mean(0)
        c = curves[name]
        lines.append(f"2. {name:8s} H(0) = {c[0]:.3f}, max {c.max():.3f} at step {int(c.argmax())}, "
                     f"H(1000) = {c[1000]:.3f}, H(3000) = {c[-1]:.3f}")
    base = np.array([run(steps=3000, seed=s, n_start=20, t_elim=None)[0] for s in range(300)]).mean(0)
    curves["fixed, no elimination"] = base
    lines.append(f"   fixed, no elimination: H(1000) = {base[1000]:.3f}, H(3000) = {base[-1]:.3f}")
    for N in (50, 100):                              # sensitivity to the number of elements, fixed regime
        Hs = np.array([run(N=N, n_start=N, steps=3000, t_elim=200, seed=s)[0] for s in range(100)]).mean(0)
        lines.append(f"   fixed, N = {N}: H(0) = {Hs[0]:.3f}, H(1000) = {Hs[1000]:.3f}, H(3000) = {Hs[-1]:.3f}, "
                     f"limit of three fit elements, log2 3 = {np.log2(3):.3f}")

    # 3. predicted against simulated change over windows of 100 steps, 60 runs per regime
    win = 100; pts = []
    for name, kw in reg.items():
        P, R = [], []
        for s in range(60):
            Hs, _, pr = run(steps=3000, seed=1000 + s, track_prediction=True, **kw)
            P.append([pr[a + 1:a + win + 1].sum() for a in range(0, 3000, win)])
            R.append([Hs[a + win] - Hs[a] for a in range(0, 3000, win)])
        P, R = np.mean(P, 0), np.mean(R, 0)
        pts += [(name, a, p, r) for a, p, r in zip(range(0, 3000, win), P, R)]
    pv = np.array([p for _, _, p, _ in pts]); rv = np.array([r for _, _, _, r in pts])
    lines.append(f"3. Windows of 100 steps, mean of 60 runs, {len(reg)} regimes ({len(pts)} windows): "
                 f"r = {np.corrcoef(pv, rv)[0, 1]:.4f}, same sign {100 * np.mean(np.sign(pv) == np.sign(rv)):.0f}%, "
                 f"least-squares slope {np.polyfit(pv, rv, 1)[0]:.3f}")

    np.savez(out / "sim_full_results.npz", steps=np.arange(3001),
             **{k.replace(" ", "_").replace(",", ""): v for k, v in curves.items()},
             pred=pv, real=rv, regime=np.array([n for n, _, _, _ in pts]))
    (out / "sim_full_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
