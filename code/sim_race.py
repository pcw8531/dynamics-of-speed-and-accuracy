"""Selection of a command as a race among the elements of the motor system (section From Averaged Information to
the Slope of the manuscript, Supplemental Material, Section S8 and Figure S5).

Each element j holds an accumulator that starts at the head start h_j and ends the selection when it reaches the
criterion theta. The accumulator of the element that the command addresses rises at the rate v, and every other
accumulator falls at the rate v, each with a variance of 2 v per unit of time. These are the mean and the variance of
the log-likelihood ratio of a signal in Gaussian noise, so the accumulators count information in natural units.

Results checked here
  1. With head start h_i the addressed element is selected after a mean time (theta - h_i) / v (Wald's identity).
  2. Another element wins with a probability of at most sum_j exp(h_j - theta). If sum_j exp(h_j) <= 1 a wrong
     selection has a probability of at most exp(-theta) whatever the shares, so one criterion holds reliability
     constant.
  3. Among all head starts within that budget, h_j = ln p_j gives the shortest mean time when commands go to
     element j with probability p_j. The mean time is then theta / v + (ln 2 / v) H, with H the averaged information
     of the shares in bits, so beta = ln 2 / v is the time in which the race gains one bit.
  4. An act of two commands selected in one race among their sequences costs theta / v + 2 beta H: the criterion is
     paid once and the head starts add.

The shares are those of the full model of practice (code/sim_full.py, fixed condition, seed 0) at seven steps.
Times are in units of 1 / v. The seed is fixed, so every run gives the same output.
Usage: python -I code/sim_race.py sim_out   (about three minutes; writes sim_out/sim_race_result.txt and
sim_out/sim_race_results.npz)
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import sim_full                                                   # noqa: E402

LN2 = np.log(2.0)
V, DT, SEED = 1.0, 0.01, 20261120                               # rate v, time step, seed
BETA_MS = 94.7 / np.log2(20)                                     # ms per bit, as in the paper
UNIT_MS = BETA_MS / LN2                                          # one unit of time (1/v) in ms


def H_bits(p):
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def race(p_cmd, starts, theta, trials, rng, chunk=5000):
    """Race among len(starts) accumulators with head starts starts and criterion theta. The addressed element is
    drawn from p_cmd.
    Returns the addressed element, the winner and the time of each trial."""
    n = len(starts)
    sd = np.sqrt(2.0 * V * DT)
    tgt_all, win_all, t_all = [], [], []
    for c0 in range(0, trials, chunk):
        m = min(chunk, trials - c0)
        tgt = rng.choice(n, size=m, p=p_cmd)
        x = np.tile(starts.astype(np.float32), (m, 1))
        drift = np.full((m, n), -V * DT, dtype=np.float32)
        drift[np.arange(m), tgt] = V * DT
        done = np.zeros(m, bool); win = np.full(m, -1); t = np.zeros(m)
        step = 0
        while not done.all():
            step += 1
            live = ~done
            x[live] += drift[live] + sd * rng.standard_normal((int(live.sum()), n), dtype=np.float32)
            hit = live & (x.max(1) >= theta)
            if hit.any():
                win[hit] = x[hit].argmax(1); t[hit] = step * DT; done |= hit
        tgt_all.append(tgt); win_all.append(win); t_all.append(t)
    return np.concatenate(tgt_all), np.concatenate(win_all), np.concatenate(t_all)


def fit(x, y):
    b, a = np.polyfit(x, y, 1)
    return a, b, float(np.corrcoef(x, y)[0, 1])


def main(out_dir):
    OUT = Path(out_dir); OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    lines = []
    def say(s=""):
        print(s); lines.append(s)

    steps = [0, 100, 200, 500, 1000, 2000, 3000]
    shares = {}
    for s in steps:
        _, k, _ = sim_full.run(N=20, n_start=20, steps=s, t_elim=200, seed=0)
        shares[s] = k / k.sum()
    say("Shares of the full model (fixed condition, seed 0)")
    for s in steps:
        say(f"  step {s:5d}  H = {H_bits(shares[s]):.3f} bits  active elements = {int((shares[s] > 0).sum())}"
            f"  largest share = {shares[s].max():.3f}")

    # ---------------- 1 to 3: one command, matched head starts, two criteria
    res = {}
    for theta in (3.0, 4.0):
        say(f"\nOne command, head starts ln p_j, criterion theta = {theta:.1f} (error bound exp(-theta) = {np.exp(-theta):.4f}), "
            f"20,000 trials per step")
        Hs, Ts, Es, I_all, T_all = [], [], [], [], []
        for s in steps:
            p = shares[s]; act = p > 0
            pa = p[act]
            tgt, win, t = race(pa, np.log(pa), theta, 20000, rng)
            ok = win == tgt
            Hs.append(H_bits(pa)); Ts.append(t.mean()); Es.append(1 - ok.mean())
            for j in range(len(pa)):
                sel = (tgt == j) & ok
                if sel.sum() >= 200:
                    I_all.append(-np.log2(pa[j])); T_all.append(t[sel].mean())
            say(f"  step {s:5d}  H = {Hs[-1]:.3f}  mean time = {Ts[-1]:.3f}  predicted theta/v + ln2 H/v = "
                f"{theta / V + LN2 * Hs[-1] / V:.3f}  error = {Es[-1]:.4f}")
        a, b, r = fit(np.array(Hs), np.array(Ts))
        say(f"  mean time on H over the 7 steps: intercept {a:.3f} (theta/v = {theta / V:.3f}), slope {b:.4f} per bit "
            f"(ln 2 / v = {LN2 / V:.4f}), r = {r:.4f}")
        a2, b2, r2 = fit(np.array(I_all), np.array(T_all))
        say(f"  time of an element on its information I_i ({len(I_all)} elements with 200 or more correct trials): "
            f"intercept {a2:.3f}, slope {b2:.4f} per bit, r = {r2:.4f}")
        say(f"  in the units of the paper (beta = {BETA_MS:.1f} ms per bit): slope {b * UNIT_MS:.1f} ms per bit, "
            f"mean time from {Ts[0] * UNIT_MS:.0f} ms at step 0 to {Ts[-1] * UNIT_MS:.0f} ms at step 3,000")
        res[theta] = dict(H=np.array(Hs), T=np.array(Ts), E=np.array(Es), I=np.array(I_all), TI=np.array(T_all))

    # ---------------- equal shares: the logarithm of the number of alternatives (Usher et al., 2002)
    say("\nEqual shares, theta = 3.0: mean time on log2 N, 20,000 trials each")
    Ns, Tn = [2, 4, 8, 16, 32], []
    for n in Ns:
        p = np.full(n, 1.0 / n)
        tgt, win, t = race(p, np.log(p), 3.0, 20000, rng)
        Tn.append(t.mean()); say(f"  N = {n:2d}  log2 N = {np.log2(n):.0f}  mean time = {t.mean():.3f}  error = {(win != tgt).mean():.4f}")
    a, b, r = fit(np.log2(Ns), np.array(Tn))
    say(f"  intercept {a:.3f}, slope {b:.4f} per bit (ln 2 / v = {LN2 / V:.4f}), r = {r:.4f}")

    # ---------------- 3: head starts that do not match the commands
    say("\nHead starts that do not match the commands, theta = 3.0, 20,000 trials each")
    p0, p3 = shares[0], shares[3000]
    a3 = p3 > 0
    cases = [("commands of step 3,000, head starts of step 3,000 (matched)", p3[a3], np.log(p3[a3])),
             ("commands of step 3,000, equal head starts", p3[a3], np.log(np.full(a3.sum(), 1.0 / a3.sum()))),
             ("commands of step 3,000, head starts of step 500", p3[a3], np.log(shares[500][a3] / shares[500][a3].sum())),
             ("commands of step 0, head starts of step 0 (matched)", p0, np.log(p0))]
    for name, pc, st in cases:
        tgt, win, t = race(pc, st, 3.0, 20000, rng)
        cross = float(-(pc * st).sum() / LN2)
        say(f"  {name}: mean time = {t.mean():.3f}, predicted theta/v + ln2 L/v = {3.0 + LN2 * cross:.3f} "
            f"with L = {cross:.3f} bits (H = {H_bits(pc):.3f}), error = {(win != tgt).mean():.4f}")

    # ---------------- 4: an act of two commands selected in one race among the 400 sequences
    say("\nAn act of two commands, one race among their sequences, theta = 3.0, 4,000 trials per step")
    H2, T2 = [], []
    for s in (0, 500, 3000):
        p = shares[s]; pa = p[p > 0]
        pp = np.outer(pa, pa).ravel()
        tgt, win, t = race(pp, np.log(pp), 3.0, 4000, rng, chunk=1000)
        H2.append(2 * H_bits(pa)); T2.append(t.mean())
        say(f"  step {s:5d}  2 H = {H2[-1]:.3f} bits  mean time = {T2[-1]:.3f}  predicted theta/v + ln2 2 H/v = "
            f"{3.0 + LN2 * H2[-1]:.3f}  error = {(win != tgt).mean():.4f}")
    Hk = np.concatenate([res[3.0]["H"], np.array(H2)]); Tk = np.concatenate([res[3.0]["T"], np.array(T2)])
    a, b, r = fit(Hk, Tk)
    say(f"  acts of one and of two commands together, mean time on the bits of the act: intercept {a:.3f}, slope {b:.4f} per bit, r = {r:.4f}")

    (OUT / "sim_race_result.txt").write_text("\n".join(lines) + "\n")
    np.savez(OUT / "sim_race_results.npz", steps=steps, H3=res[3.0]["H"], T3=res[3.0]["T"], E3=res[3.0]["E"],
             I3=res[3.0]["I"], TI3=res[3.0]["TI"], H4=res[4.0]["H"], T4=res[4.0]["T"], E4=res[4.0]["E"],
             Ns=Ns, Tn=Tn, H2=H2, T2=T2)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
