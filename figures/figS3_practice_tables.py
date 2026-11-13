"""Figure S3 of the Supplemental Material: the two classical practice tables of the slope test.

A, the classical law at each stage of Leonard's (1954) practice as tabled by Mowbray and Rhoades (1959, Table I),
three and six choices in May, September and December, with the slope of each month in ms per bit. B, the same for the
participant of Mowbray and Rhoades (1959, Table II), two and four choices, each trial about 500 reactions, with single
trials (thin lines) and the mean of trials 1 to 3 and of trials 13 to 15. Both pairs of set sizes differ by one bit, so
each slope is a plain difference in milliseconds. C, time saved by practice, from the first stage to the last, against
the information of the choice (Equation 17). The reading of these slopes as averaged information is Figure 10 of the
manuscript (figures/fig10_slope_reading.py).
Inputs: data/transcribed/mowbray_rhoades1959_table2.csv, data/transcribed/mowbray_rhoades1959_tables1_3.csv.
Run from the repository root: python figures/figS3_practice_tables.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import Sheet, INK, GREY, TONE, CLEAR, tb_axes, save_panels  # noqa: E402

# ---------------------------------------------------------------- data
mr = pd.read_csv(ROOT / "data" / "transcribed" / "mowbray_rhoades1959_table2.csv")
t13 = pd.read_csv(ROOT / "data" / "transcribed" / "mowbray_rhoades1959_tables1_3.csv")
MONTHS = ["May", "September", "December"]
lp = t13[t13.table == "I"].pivot(index="period_or_trial", columns="choices", values="mean_rt_ms").loc[MONTHS]
mp = mr.pivot(index="trial", columns="choices", values="mean_rt_ms")
leo_b = (lp[6] - lp[3]).values                                      # ms per bit, log2 6 - log2 3 = 1
mr_b = (mp[4] - mp[2])                                              # ms per bit, log2 4 - log2 2 = 1
BLOCKS = {"trials 1 to 3": [1, 2, 3], "trials 13 to 15": [13, 14, 15]}
mr_blk = {k: mp.loc[v].mean() for k, v in BLOCKS.items()}           # block means of each set size
mr_blk_b = {k: float(v[4] - v[2]) for k, v in mr_blk.items()}
saved = {"Leonard": (lp[3].iloc[0] - lp[3].iloc[-1], lp[6].iloc[0] - lp[6].iloc[-1]),
         "Mowbray": (mp[2].iloc[0] - mp[2].iloc[-1], mp[4].iloc[0] - mp[4].iloc[-1])}
assert list(leo_b) == [62, 51, 45] and abs(mr_blk_b["trials 1 to 3"] - 18.67) < 0.01 and mr_blk_b["trials 13 to 15"] == 2.0
assert all(lg > sm for sm, lg in saved.values())

# ---------------------------------------------------------------- layout
S = Sheet(6.5, 2.0)
Y1, HH, WW = 0.45, 1.27, 1.50
axA, axB, axC = S.ax(0.52, Y1, WW, HH), S.ax(2.66, Y1, WW, HH), S.ax(4.75, Y1, WW, HH)
FIRST, MID, LAST = TONE["small-world"], TONE["scale-free"], INK
EQ = dict(fontsize=6.4)

# A, Leonard
ax = axA
x3, x6 = np.log2(3), np.log2(6)
for m, c, b_, dy in zip(MONTHS, (FIRST, MID, LAST), leo_b, (6, -4, 0)):
    ax.plot([x3, x6], [lp.loc[m, 3], lp.loc[m, 6]], "-o", color=c, ms=3.0, lw=1.1, mec="white", mew=0.4)
    ax.text(x6 + 0.06, lp.loc[m, 6] + dy, f"{m[:3]}, {b_:.0f}", va="center", ha="left", color=c, **EQ, **CLEAR)
ax.set_xlim(1.35, 3.15); ax.set_xticks([x3, x6]); ax.set_xticklabels(["3", "6"])
tb_axes(ax, yticks=[250, 290, 330], ylim=(235, 335))
ax.set_xlabel("choices $n$, on $\\log_2 n$", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)
ax.text(1.40, 328, "slope, ms per bit", color=GREY, **EQ, **CLEAR)

# B, Mowbray and Rhoades
ax = axB
for t in mp.index:
    c = FIRST if t <= 3 else LAST
    ax.plot([1, 2], [mp.loc[t, 2], mp.loc[t, 4]], color=c, lw=0.5, alpha=0.6)
for (k, v), c in zip(mr_blk.items(), (FIRST, LAST)):
    ax.plot([1, 2], [v[2], v[4]], "-o", color=c, ms=3.0, lw=1.6, mec="white", mew=0.4)
ax.text(2.06, mr_blk["trials 1 to 3"][4] + 3, f"trials 1 to 3, {mr_blk_b['trials 1 to 3']:.0f}", va="center", ha="left", color=FIRST, **EQ)
ax.text(2.06, mr_blk["trials 13 to 15"][4] - 6, f"13 to 15, {mr_blk_b['trials 13 to 15']:.0f}", va="center", ha="left", color=LAST, **EQ)
ax.set_xlim(0.85, 2.95); ax.set_xticks([1, 2]); ax.set_xticklabels(["2", "4"])
tb_axes(ax, yticks=[210, 260, 310], ylim=(205, 315))
ax.set_xlabel("choices $n$, on $\\log_2 n$", labelpad=1); ax.set_ylabel("reaction time, ms", labelpad=1.5)
ax.text(0.90, 308, "thin lines, single trials", color=GREY, **EQ, **CLEAR)

# C, time saved against information (Equation 17)
ax = axC
for (lab, (sm, lg)), xs, c, mk in zip(saved.items(), ((x3, x6), (1.0, 2.0)), (MID, LAST), ("o", "s")):
    ax.plot(xs, [sm, lg], "-", color=c, lw=1.1, zorder=2)
    ax.plot(xs, [sm, lg], mk, color=c, ms=3.4, mec="white", mew=0.4, zorder=3)
    for x_, y_ in zip(xs, (sm, lg)):
        ax.text(x_, y_ + 4, f"{y_:.0f}", ha="center", va="bottom", fontsize=6.2, color=c)
ax.text(1.80, 6, f"Leonard, {saved['Leonard'][1] - saved['Leonard'][0]:.0f} ms per bit", color=MID, **EQ)
ax.text(1.02, 60, f"Mowbray and Rhoades, {saved['Mowbray'][1] - saved['Mowbray'][0]:.0f} ms per bit", color=LAST, **EQ)
tb_axes(ax, xticks=[1, 2, 3], xlim=(0.8, 3.05), yticks=[0, 50, 100], ylim=(0, 108))
ax.set_xlabel(r"information of the choice, $\log_2 n$, bits", labelpad=1); ax.set_ylabel("time saved, ms", labelpad=1.5)

S.head_over(axA, "A", "Leonard (1954)")
S.head_over(axB, "B", "Mowbray and Rhoades (1959)")
S.head_over(axC, "C", "Time saved grows with difficulty")
out = ROOT / "figures" / "figS3_practice_tables.png"
S.fig.savefig(out)
(ROOT / "figures" / "panels").mkdir(exist_ok=True)
save_panels(S, {"A": [axA], "B": [axB], "C": [axC]}, str(ROOT / "figures" / "panels" / "figS3"))
print("Leonard slopes", leo_b.tolist(), "| Mowbray and Rhoades block slopes", {k: round(v, 2) for k, v in mr_blk_b.items()},
      "| slopes by trial", mr_b.values.tolist())
print("time saved (smaller, larger):", {k: tuple(float(x) for x in v) for k, v in saved.items()})
print("written", out)
