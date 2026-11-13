"""The classical laws as static cases of the dynamic law, by simulation (Supplemental Section S4.4).

One simulated performer does three tasks twice, with the coupling of the full model at step 0 (even coupling over 20
elements, H = log2 20) and at step 1,000 (one run of Figure 5, chosen as the run whose H at step 1,000 is nearest the
mean of 300 runs). Nothing in the tasks is fitted at step 1,000, so the times there are predictions.

Commands.  A command addresses one element, drawn by its share of coupling p_i. Commands are coded three at a time by a
Huffman code over all triples, so a command takes L bits on average with H <= L < H + 1/3 (Equation 9). Each
command bit takes beta ms.

Hick.     A choice among n equally likely lights, n = 1 to 8, is made by halvings. The lights of three trials are coded
          together by a Huffman code over all triples, and each binary digit of the codeword is one halving, so a trial
          needs log2 n halvings on average to within 1/3. Each halving issues one command.
Hyman.    The same, with the probabilities of the 24 conditions of Hyman (1953), Tables 1 and 2. In Experiment III the
          lights are coded by their probabilities after the previous light, which have the same shape after every light.
          Only condition means are formed.
Fitts.    An aimed movement is an initial adjustment of a0 = 180 ms, the time at ID 1 in Fitts (1954, Table 1), which
          leaves log2(D/E0) = (a0 - a) / b = 1.766 bits fewer to remove, then corrections of one bit each, one command
          each, until the error is within W/2. For ID - 1.766 = s + f, a trial makes s corrections and one more with
          probability f; below 1.766 bits it makes none.

Calibration at step 0.  Hick and Hyman: a and beta by least squares of the observed times of participant G. C.
(Hyman, 1953, Figure 1) on the expected command bits per trial, Experiment I for Hick and all 24 conditions for Hyman.
Fitts: a0 = 180 ms, beta = b / L(0) with b = 94.7 ms per bit, the 1-oz line of Table 1.

Each condition is run for T trials and the whole set is repeated R times with new draws. The slope of each task is the
least-squares slope of its condition means on the information of the task (log2 n, the information of the condition,
ID from 2 to 7 bits), and the ratio of the slope at step 1,000 to the slope at step 0 is set against H(1000)/H(0).

Usage: python -I code/sim_static_cases.py <output folder>
"""
import heapq
import itertools
import re
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sim_full as sf  # noqa: E402

R_BLOCK, T, R, SEED_RUNS = 3, 3000, 200, 300
A0_FITTS, B_FITTS, A_FITTS = 180.0, 94.7, 12.8


# ---------------------------------------------------------------- coding
def huffman_lengths(p):
    """Codeword lengths of a binary Huffman code for probabilities p (all > 0). One symbol needs no digit."""
    p = np.asarray(p, float)
    if len(p) == 1:
        return np.zeros(1, int)
    heap = [(float(q), i, (i,)) for i, q in enumerate(p)]
    heapq.heapify(heap); L = np.zeros(len(p), int); tie = len(p)
    while len(heap) > 1:
        q1, _, s1 = heapq.heappop(heap); q2, _, s2 = heapq.heappop(heap)
        for i in s1 + s2:
            L[i] += 1
        heapq.heappush(heap, (q1 + q2, tie, s1 + s2)); tie += 1
    return L


class BlockCode:
    """Huffman code over blocks of r symbols drawn independently from p."""

    def __init__(self, p, r=R_BLOCK):
        p = np.asarray(p, float); p = p[p > 0] / p[p > 0].sum()
        self.p, self.r = p, r
        P = np.array([np.prod(c) for c in itertools.product(p, repeat=r)])
        self.P, self.len = P / P.sum(), huffman_lengths(P / P.sum())
        self.cum = np.cumsum(self.P)
        self.len1 = huffman_lengths(p)                       # for a remainder shorter than a block
        self.H = float(-(p * np.log2(p)).sum())
        self.L = float((self.P * self.len).sum()) / r        # digits per symbol, H <= L < H + 1/r

    def digits(self, n, rng):
        """Total digits for n symbols drawn from p, coded r at a time."""
        nb, rem = divmod(int(n), self.r)
        u = rng.random(nb)
        tot = int(self.len[np.minimum(np.searchsorted(self.cum, u, side="right"), len(self.P) - 1)].sum())
        if rem:
            idx = np.minimum(np.searchsorted(np.cumsum(self.p), rng.random(rem), side="right"), len(self.p) - 1)
            tot += int(self.len1[idx].sum())
        return tot


# ---------------------------------------------------------------- tasks
def parse_probs(text, lights):
    """Probabilities of one condition from Hyman's tables, e.g. 'after a: a 7/10, each other 1/10'."""
    text = text.split(":", 1)[1] if ":" in text else text
    if "," not in text and " " in text.strip() and not re.search(r"[a-z]", text):
        return [Fraction(x) for x in text.split()]
    named, other = [], None
    for part in text.split(","):
        part = part.strip()
        if part.startswith("each other"):
            other = Fraction(part.split()[-1])
        elif re.match(r"^[a-h] ", part):
            named.append(Fraction(part.split()[-1]))
        else:
            return [Fraction(x) for x in part.split()]
    probs = named + ([other] * (lights - len(named)) if other is not None else [])
    assert len(probs) == lights and sum(probs) == 1, (text, probs)
    return [q for q in probs if q > 0]


def hyman_conditions():
    """The 24 conditions of participant G. C.: experiment, condition, probabilities, information, observed RT."""
    tab = pd.read_csv(ROOT / "data" / "transcribed" / "hyman1953_tables1_2.csv")
    e1 = pd.read_csv(ROOT / "data" / "read_from_figures" / "hyman1953_fig1_GC_experiment1.csv")
    e23 = pd.read_csv(ROOT / "data" / "read_from_figures" / "hyman1953_fig1_GC_experiments2_3.csv")
    rows = [("I", int(n), [1.0 / n] * int(n), float(np.log2(n)), float(rt)) for n, rt in zip(e1.alternatives, e1.rt_ms)]
    for _, t in tab.iterrows():
        rt = float(e23[(e23.experiment == t.experiment) & (e23.condition == t.condition)].rt_ms.iloc[0])
        q = [float(x) for x in parse_probs(t.probabilities, int(t.lights))]
        Hq = float(-(np.array(q) * np.log2(q)).sum())
        assert abs(Hq - t.information_bits) < 0.011, (t.experiment, t.condition, Hq)
        rows.append((t.experiment, int(t.condition), q, float(t.information_bits), rt))
    return rows


def condition_mean(halvings_per_trial_code, command_code, a, beta, rng, fitts_bits=None):
    """Mean time of T trials. Hick and Hyman: halvings from the block code of the lights. Fitts: corrections."""
    if fitts_bits is None:
        n_cmd = halvings_per_trial_code.digits(T, rng)
    else:
        x = max(0.0, fitts_bits); s, f = int(np.floor(x)), x - np.floor(x)
        n_cmd = T * s + int((rng.random(T) < f).sum())
    return a + beta * command_code.digits(n_cmd, rng) / T


def coupling_at(step):
    """Shares of coupling of the chosen run at a step: step 0 is the even ring."""
    finals = np.array([sf.run(N=20, n_start=20, steps=1000, lam=0.5, t_elim=200, seed=s)[0][-1]
                       for s in range(SEED_RUNS)])
    seed = int(np.argmin(np.abs(finals - finals.mean())))
    _, k, _ = sf.run(N=20, n_start=20, steps=step, lam=0.5, t_elim=200, seed=seed)
    return k / k.sum(), seed, float(finals.mean())


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    p1000, seed, H_mean = coupling_at(1000)
    p0 = np.full(20, 1 / 20)
    code = {0: BlockCode(p0), 1000: BlockCode(p1000)}
    lines = [f"Coupling at step 1,000: run {seed} of {SEED_RUNS}, H = {code[1000].H:.4f} bits "
             f"(mean of the runs {H_mean:.4f}); step 0, H = {code[0].H:.4f}",
             f"Command bits, blocks of {R_BLOCK}: L(0) = {code[0].L:.4f}, L(1000) = {code[1000].L:.4f}; "
             f"H(1000)/H(0) = {code[1000].H / code[0].H:.4f}, L(1000)/L(0) = {code[1000].L / code[0].L:.4f}"]
    for st in (0, 1000):
        assert code[st].H <= code[st].L < code[st].H + 1 / R_BLOCK

    hy = hyman_conditions()
    light = [BlockCode(q) for _, _, q, _, _ in hy]
    for c, (_, _, _, h, _) in zip(light, hy):
        assert c.H - 1e-12 <= c.L < c.H + 1 / R_BLOCK
    x_hy = np.array([h for _, _, _, h, _ in hy]); rt_hy = np.array([r for *_, r in hy])
    ex1 = np.array([e == "I" for e, *_ in hy])
    halv = np.array([c.L for c in light])                   # expected halvings per trial

    # calibration at step 0 on expected command bits per trial
    b_hick, a_hick = np.polyfit(halv[ex1] * code[0].L, rt_hy[ex1], 1)
    b_hym, a_hym = np.polyfit(halv * code[0].L, rt_hy, 1)
    fitts = pd.read_csv(ROOT / "data" / "transcribed" / "fitts1954_table1.csv")
    ID = np.arange(1, 8, dtype=float); c_fitts = (A0_FITTS - A_FITTS) / B_FITTS
    beta = {"Hick": b_hick, "Hyman": b_hym, "Fitts": B_FITTS / code[0].L}
    lines.append(f"Calibration at step 0: Hick a = {a_hick:.1f} ms, beta = {b_hick:.2f} ms per command bit; "
                 f"Hyman a = {a_hym:.1f}, beta = {b_hym:.2f}; Fitts a0 = {A0_FITTS:.0f}, log2(D/E0) = {c_fitts:.3f}, "
                 f"beta = {beta['Fitts']:.2f}")

    rng = np.random.default_rng(20261107)
    means = {(t, st): [] for t in ("Hick", "Hyman", "Fitts") for st in (0, 1000)}
    for _ in range(R):
        for st in (0, 1000):
            means[("Hick", st)].append([condition_mean(c, code[st], a_hick, b_hick, rng)
                                        for c, e in zip(light, ex1) if e])
            means[("Hyman", st)].append([condition_mean(c, code[st], a_hym, b_hym, rng) for c in light])
            means[("Fitts", st)].append([condition_mean(None, code[st], A0_FITTS, beta["Fitts"], rng,
                                                        fitts_bits=i - c_fitts) for i in ID])
    means = {k: np.array(v) for k, v in means.items()}
    xs = {"Hick": x_hy[ex1], "Hyman": x_hy, "Fitts": ID}
    use = {"Hick": slice(None), "Hyman": slice(None), "Fitts": ID >= 2}
    slope = {k: np.array([np.polyfit(xs[k[0]][use[k[0]]], m[use[k[0]]], 1)[0] for m in v]) for k, v in means.items()}
    ratio = {t: slope[(t, 1000)] / slope[(t, 0)] for t in ("Hick", "Hyman", "Fitts")}
    for t in ("Hick", "Hyman", "Fitts"):
        lo, hi = np.percentile(ratio[t], [2.5, 97.5])
        lines.append(f"{t:5s} slope at step 0 {slope[(t, 0)].mean():6.1f}, at step 1,000 {slope[(t, 1000)].mean():6.1f} "
                     f"ms per bit; ratio {ratio[t].mean():.4f} (95% of repetitions {lo:.4f} to {hi:.4f})")
    obs = {"Hick": rt_hy[ex1], "Fitts": 1000 * fitts.groupby("Id_bits").t_1oz_s.mean().reindex(ID.astype(int)).values}
    lines.append(f"T = {T} trials per condition, R = {R} repetitions, blocks of {R_BLOCK}")

    np.savez(out / "sim_static_cases.npz", x_hick=xs["Hick"], x_hyman=xs["Hyman"], x_fitts=ID,
             obs_hick=obs["Hick"], obs_fitts=obs["Fitts"], fitts_id=fitts.Id_bits.values,
             fitts_mt=1000 * fitts.t_1oz_s.values, experiment=np.array([e for e, *_ in hy]),
             H0=code[0].H, H1000=code[1000].H, L0=code[0].L, L1000=code[1000].L, p1000=p1000, seed=seed,
             **{f"mean_{t}_{st}": v for (t, st), v in means.items()},
             **{f"ratio_{t}": v for t, v in ratio.items()})
    (out / "sim_static_cases_result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim_out")
