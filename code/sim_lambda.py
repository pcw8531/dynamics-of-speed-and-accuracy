"""How strongly practice is directed: the fixed condition of the full model (sim_full.py) at four values of lambda.

lambda = 0 is practice without direction (every element has fitness 1), lambda > 0 gives the three targeted elements
fitness 1 + lambda (Equation 13). 300 runs each, seeds 0 to 299, to step 1,000.
Usage: python -I code/sim_lambda.py   (writes sim_out/sim_lambda_result.txt)
"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import sim_full as sf                                                 # noqa: E402

lines = ["Fixed condition of the full model, N = 20, elimination from step 200, 300 runs, mean H in bits (and the slope beta H, ms per bit)"]
beta = 94.7 / np.log2(20)
for lam in (0.0, 0.25, 0.5, 1.0):
    H = np.array([sf.run(N=20, n_start=20, steps=1000, lam=lam, t_elim=200, seed=s)[0] for s in range(300)]).mean(0)
    lines.append(f"lambda = {lam:<4}: H(0) = {H[0]:.3f}, H(200) = {H[200]:.3f}, H(500) = {H[500]:.3f}, H(1000) = {H[1000]:.3f}"
                 f"  (slope {beta * H[0]:.1f} -> {beta * H[1000]:.1f})")
(ROOT / "sim_out" / "sim_lambda_result.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
