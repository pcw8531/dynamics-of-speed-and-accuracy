"""Check that preferential attachment on a fixed set of elements settles at the Polya-urn limit (Supplemental Material, Section S3.2).

Independent draws: each link end is drawn with probability k_i / K and adds 1 to k_i. The shares converge
to a Dirichlet distribution with parameters k_i(0), whose expected entropy is
E[H_nat] = psi(A + 1) - sum_i (a_i / A) psi(a_i + 1), A = sum_i a_i.
Distinct ends: the rule of sim_practice.py, where the two ends of a link fall on different elements.
"""
import numpy as np
from scipy.special import digamma

N, A0 = 20, 2.0
RUNS, ENDS = 400, 8000


def predicted_bits(n=N, a0=A0):
    return (digamma(n * a0 + 1) - digamma(a0 + 1)) / np.log(2)


def final_entropy(rng, distinct):
    k = np.full(N, A0)
    steps = ENDS // 2 if distinct else ENDS
    for _ in range(steps):
        c = np.cumsum(k)
        i = np.searchsorted(c, rng.random() * c[-1], side="right")
        if distinct:
            w = k.copy(); w[i] = 0.0
            c2 = np.cumsum(w)
            j = np.searchsorted(c2, rng.random() * c2[-1], side="right")
            k[i] += 1; k[j] += 1
        else:
            k[i] += 1
    p = k / k.sum()
    return float(-(p * np.log2(p)).sum())


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    print(f"Dirichlet prediction: {predicted_bits():.4f} bits")
    for distinct in (False, True):
        H = np.array([final_entropy(rng, distinct) for _ in range(RUNS)])
        label = "distinct ends    " if distinct else "independent draws"
        print(f"{label}: {H.mean():.4f} +/- {H.std(ddof=1) / np.sqrt(RUNS):.4f} bits ({RUNS} runs, {ENDS} link ends)")
