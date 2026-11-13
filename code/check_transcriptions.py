"""Internal consistency checks on transcribed tables and on values read from figures.
Prints every check and exits non-zero on failure."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent.parent / "data" / "transcribed"
RF = T.parent / "read_from_figures"
fail = 0

f = pd.read_csv(T / "fitts1954_table1.csv")
id_calc = np.log2(2 * f.A_in / f.W_in)
ok = np.allclose(id_calc, f.Id_bits)
print(f"Fitts Table 1: Id = log2(2A/W) in all {len(f)} rows: {ok}")
fail += not ok
for s in ("1oz", "1lb"):
    ip = f.Id_bits / f[f"t_{s}_s"]
    dev = (ip - f[f"Ip_{s}_bits_per_s"]).abs()
    ok = (dev <= 0.03).all()
    print(f"Fitts Table 1, {s}: Id/t matches printed Ip within 0.03 bits/s: {ok} (largest gap {dev.max():.3f})")
    fail += not ok

b, a = np.polyfit(f.Id_bits, 1000 * f.t_1oz_s, 1)
r = np.corrcoef(f.Id_bits, f.t_1oz_s)[0, 1]
ok = round(a, 1) == 12.8 and round(b, 1) == 94.7 and round(r, 4) == 0.9831
print(f"Fitts Table 1, 1oz: least-squares line MT = {a:.1f} + {b:.1f} ID ms, r = {r:.4f}; matches MacKenzie (1992): {ok}")
fail += not ok

def rank_check(df, ip, rank):
    expected = df[ip].rank(ascending=False, method="average")
    return bool((expected == df[rank]).all())

print(f"Fitts Table 1: printed ranks follow Ip (1oz, 1lb): {rank_check(f, 'Ip_1oz_bits_per_s', 'rank_1oz')}, {rank_check(f, 'Ip_1lb_bits_per_s', 'rank_1lb')}")
for name in ("fitts1954_table2.csv", "fitts1954_table3.csv"):
    g = pd.read_csv(T / name)
    ok1 = np.allclose(np.log2(2 * g.A_in / g.W_in), g.Id_bits)
    gap = ((g.Id_bits / g.t_s) - g.Ip_bits_per_s).abs().max()
    ok2 = gap <= 0.03
    ok3 = rank_check(g, "Ip_bits_per_s", "rank")
    print(f"{name}: Id = log2(2A/W) {ok1}; Id/t matches Ip within 0.03 {ok2} (largest gap {gap:.3f}); ranks follow Ip {ok3}")
    fail += not (ok1 and ok2 and ok3)

l = pd.read_csv(T / "mowbray_rhoades1959_tables1_3.csv")
t1 = l[l.table == "I"].pivot(index="period_or_trial", columns="choices", values="mean_rt_ms")
diff = (t1[6] - t1[3]).reindex(["May", "September", "December"]).tolist()
ok = diff == [62, 51, 45]
print(f"Mowbray and Rhoades Table I: six minus three choices = {diff}, matches printed 62, 51, 45: {ok}")
fail += not ok

m = pd.read_csv(T / "mowbray_rhoades1959_table2.csv")
ok = len(m) == 12 and set(m.trial) == {1, 2, 3, 13, 14, 15} and set(m.choices) == {2, 4}
print(f"Mowbray and Rhoades Table II: 12 cells, trials 1-3 and 13-15, two and four choices: {ok}")
fail += not ok

w = pd.read_csv(T / "woodworth1899_tables_I-IV_right_hand.csv")
ok = len(w) == 80 and w.average_error_mm.isna().sum() == 2
print(f"Woodworth Tables I-IV: 80 cells, 2 blank as in the original: {ok}")
fail += not ok

h = pd.read_csv(RF / "hyman1953_fig1_GC_experiment1.csv")
gapAB = (h.rt_ms_reading_A - h.rt_ms_reading_B).abs().max()
gapx = (h.bits_read_from_figure - h.information_bits).abs().max()
r = np.corrcoef(h.information_bits, h.rt_ms)[0, 1]
ok = gapAB <= 2 and gapx <= 0.02 and abs(r - 0.991) <= 0.005
print(f"Hyman Figure 1, G.C., Exp. I: readings A and B within {gapAB:.1f} ms; bits read within {gapx:.3f} of text; r = {r:.3f} vs .991 printed: {ok}")
fail += not ok

h23 = pd.read_csv(RF / "hyman1953_fig1_GC_experiments2_3.csv")
tb = pd.read_csv(T / "hyman1953_tables1_2.csv")
gapAB = (h23.rt_ms_reading_A - h23.rt_ms_reading_B).abs().max()
gapx = (h23.bits_read_from_figure - h23.information_bits).abs().max()
same = h23.merge(tb, on=["experiment", "condition"], suffixes=("", "_table"))
ok_bits = len(same) == 16 and (same.information_bits == same.information_bits_table).all()
rs = {e: np.corrcoef(g.information_bits, g.rt_ms)[0, 1] for e, g in h23.groupby("experiment")}
allp = pd.concat([h[["information_bits", "rt_ms"]], h23[["information_bits", "rt_ms"]]])
rj = np.corrcoef(allp.information_bits, allp.rt_ms)[0, 1]
ok = gapAB <= 3.5 and gapx <= 0.06 and ok_bits and abs(rs["II"] - 0.989) <= 0.005 and abs(rs["III"] - 0.972) <= 0.005 \
    and abs(rj - 0.985) <= 0.003
print(f"Hyman Figure 1, G.C., Exp. II and III: readings A and B within {gapAB:.1f} ms; bits read within {gapx:.3f} of Tables 1"
      f" and 2; r = {rs['II']:.3f} and {rs['III']:.3f} vs .989 and .972 printed; all 24, r = {rj:.3f} vs .985 printed: {ok}")
fail += not ok

from fractions import Fraction
ok = True
for _, row in tb.iterrows():
    fr = [float(Fraction(x)) for x in __import__("re").findall(r"\d+/\d+", row.probabilities)]
    if row.experiment == "II":
        q = fr
    elif row.probabilities.startswith("after a: a 0"):
        q = [fr[0]] * (row.lights - 1)
    elif "each other" in row.probabilities:
        q = [fr[0]] + [fr[1]] * (row.lights - 1)
    else:
        q = fr
    Hq = -sum(x * np.log2(x) for x in q)
    ok &= abs(sum(q) - 1) < 1e-9 and abs(Hq - row.information_bits) < 0.011
print(f"Hyman Tables 1 and 2: 16 conditions, probabilities sum to 1 and give the printed bits within 0.01: {ok}")
fail += not ok

c = pd.read_csv(T / "hyman1953_p194_components.csv").set_index("component")
wt = {"frequent light": 13 / 16, "each rare light": 3 / 16}
obs = sum(w * c.loc[k, "rt_observed_ms"] for k, w in wt.items())
pre = sum(w * c.loc[k, "rt_predicted_ms"] for k, w in wt.items())
line = [212 + 153 * c.loc[k, "information_bits"] for k in wt]
# weighting the printed component means by 13/16 and 3/16 reproduces the printed condition means within 3 ms
ok = abs(obs - 361) < 3.5 and abs(pre - 363) < 1.5 and all(abs(l_ - c.loc[k, "rt_predicted_ms"]) < 3 for l_, k in zip(line, wt))
print(f"Hyman p. 194, Exp. II Cond. 3: weighted means {obs:.1f} observed and {pre:.1f} predicted vs 361 and 363 printed; "
      f"predictions on 212 + 153 H: {ok}")
fail += not ok

z = pd.read_csv(RF / "schmidt1979_fig8_fig9_distance_We.csv")
gap = (z.We_mm_fig8 - z.We_mm_fig9).abs().max()
r = np.corrcoef(z.A_over_MT_fig9_cm_s, z.We_mm_fig9)[0, 1]
ok = len(z) == 9 and gap <= 0.2 and 0.965 <= r < 0.985
print(f"Schmidt et al. Figures 8 and 9: 9 means; two readings of We within {gap:.2f} mm; r = {r:.3f} vs .97 printed: {ok}")
fail += not ok

sys.exit(1 if fail else 0)
