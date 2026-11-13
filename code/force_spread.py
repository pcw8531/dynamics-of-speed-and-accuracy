"""The impulse-variability equation as a static case of the dynamic law, and its slope at two steps of practice.

Equation 22 of the main text (Supplemental Section S4.5) gives the relative spread of the force of an act shared among
the elements in proportion to their shares p_i, each element's force varying with a standard deviation kappa times its
mean and r_ij the correlation of the errors of elements i and j:

    sigma_F / F = kappa [sum_i p_i^2 + sum_{i != j} p_i p_j r_ij]^(1/2) = kappa [S + (1 - S) rbar]^(1/2),

with S = sum_i p_i^2 and rbar = sum_{i != j} p_i p_j r_ij / (1 - S), the mean correlation of the errors weighted by the
shares (an identity, since sum_{i != j} p_i p_j = 1 - S). The slope b(t) of Equation 3 follows sigma_F / F, so
b(t) / b(0) = sqrt[(S_t + (1 - S_t) rbar_t) / (S_0 + (1 - S_0) rbar_0)].

The shares are those of the simulated performer of Section S4.4 (sim_static_cases.py): even over 20 elements at step 0,
and at step 1,000 the run of the full model whose H lies nearest the mean of 300 runs (sim_out/sim_static_cases.npz,
p1000). At step 0 the errors are taken as independent (rbar = 0), and the line of step 0 is the least-squares line
through the nine effective widths read from Schmidt et al. (1979, Figure 9), as in Figure 1D and Table 2. The script
gives, at step 1,000, the slope with independent errors, the rbar that keeps the slope of step 0, the rbar that lowers it
by the common factor of the other tasks (the ratio of command bits of the simulated performer, 0.697), and the lowest
rbar that a valid covariance allows at these shares, -S / (1 - S), where the summed force would not vary at all.

Run from the repository root: python code/force_spread.py sim_out
Writes sim_out/force_spread_result.txt.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
N = 20


def spread(S, rbar):
    """sigma_F / F in units of kappa."""
    return np.sqrt(S + (1.0 - S) * rbar)


def rbar_for(q, S, S0=1.0 / N):
    """The share-weighted correlation at which the slope is q times that of even shares with independent errors."""
    return (q * q * S0 - S) / (1.0 - S)


def main(out_dir):
    out = Path(out_dir)
    d = np.load(out / "sim_static_cases.npz")
    p = d["p1000"]; S0, S1 = 1.0 / N, float((p ** 2).sum())
    H1 = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    q = float(d["L1000"] / d["L0"])                                  # the common factor of the three time slopes
    z = pd.read_csv(ROOT / "data" / "read_from_figures" / "schmidt1979_fig8_fig9_distance_We.csv")
    b0, a0 = np.polyfit(z.A_over_MT_fig9_cm_s, z.We_mm_fig9, 1)
    r0 = np.corrcoef(z.A_over_MT_fig9_cm_s, z.We_mm_fig9)[0, 1]
    lines = [
        f"Line through the nine values (Schmidt et al., 1979, Figure 9): W_e = {a0:.3f} + {b0:.4f} D/MT, r = {r0:.4f}",
        f"Step 0: even shares over {N} elements, sum p^2 = {S0:.4f}, H_2 = {-np.log2(S0):.3f} bits, "
        f"sigma_F/F = {spread(S0, 0):.4f} kappa (independent errors)",
        f"Step 1,000: run {int(d['seed'])}, H = {H1:.4f} bits, sum p^2 = {S1:.4f}, H_2 = {-np.log2(S1):.3f} bits",
        f"  independent errors (rbar = 0): sigma_F/F = {spread(S1, 0):.4f} kappa, slope ratio "
        f"{spread(S1, 0) / spread(S0, 0):.3f}, slope {b0 * spread(S1, 0) / spread(S0, 0):.4f}",
        f"  slope of step 0 kept at rbar = {rbar_for(1.0, S1):.4f}",
        f"  slope lowered by the common factor {q:.4f} at rbar = {rbar_for(q, S1):.4f}, slope {b0 * q:.4f}, "
        f"variance of the summed force {q * q * S0 / S1:.3f} of that with independent errors",
        f"  lowest rbar of a valid covariance at these shares: {-S1 / (1 - S1):.4f}",
    ]
    (out / "force_spread_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
