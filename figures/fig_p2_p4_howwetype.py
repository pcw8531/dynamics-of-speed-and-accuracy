"""Per-typist values from the Prediction 2 and Prediction 4 tests in How We Type, as the plans require.
Reads analysis/p4_per_typist.csv. Run from the repository root."""
import logging
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
ROOT = Path(__file__).resolve().parent.parent
p = pd.read_csv(ROOT / "analysis" / "p4_per_typist.csv")
plt.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "axes.linewidth": 0.7})
fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.3)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.2, top=0.88, wspace=0.42)
ceil = np.log2(50)
for ax, col, lab, t in zip(axs, ["slope_ms_per_bit", "median_iki_ms", "corrections_per_char"],
                           ["slope (ms per bit)", "median inter-key interval (ms)", "corrections per character"],
                           ["A  Slope", "B  Speed", "C  Corrections"]):
    for g, mk, fc in (("faster", "o", "black"), ("slower", "s", "white")):
        d = p[p.group == g]
        ax.plot(d.H_bits, d[col], mk, ms=3.6, mfc=fc, mec="black", mew=0.7, label=g)
    ax.axvline(ceil, color="#6e6e6e", lw=0.7, ls=":")
    ax.set_xlabel("H (bits)"); ax.set_ylabel(lab); ax.set_title(t, loc="left", fontsize=9, fontweight="bold")
    ax.set_xlim(5.605, 5.648); ax.set_xticks([5.61, 5.63, 5.644]); ax.set_xticklabels(["5.61", "5.63", "log$_2$50"])
axs[0].legend(frameon=False, fontsize=7, loc="lower left", handletextpad=0.2)
fig.savefig(ROOT / "figures" / "fig_p2_p4_howwetype.png"); print("written")
