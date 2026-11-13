"""Session mean H and ball speed across the 15 practice sessions, each participant and side (analysis/throwing_per_session.csv)."""
import logging
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
ROOT = Path(__file__).resolve().parent.parent
s = pd.read_csv(ROOT / "analysis" / "throwing_per_session.csv")
plt.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "axes.linewidth": 0.7})
fig, axs = plt.subplots(2, 2, figsize=(6.4, 4.4), sharex=True); fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.11, hspace=0.3, wspace=0.28)
greys = ["#000000", "#3a3a3a", "#5a5a5a", "#787878", "#959595", "#b0b0b0", "#c8c8c8"]
for j, side in enumerate(["DP", "NDP"]):
    for i, (col, lab) in enumerate([("H", "session mean H (bits)"), ("speed", "ball speed")]):
        ax = axs[i, j]
        for k, subj in enumerate(range(1, 8)):
            g = s[(s.subject == subj) & (s.side == side)]
            ax.plot(g.session_no, g[col], "-o", ms=2.2, lw=0.8, color=greys[k], label=f"S{subj}")
        ax.set_ylabel(lab)
        if i == 1: ax.set_xlabel("session")
        ax.set_title(f"{'ABCD'[2 * i + j]}  {'dominant' if side == 'DP' else 'non-dominant'} arm", loc="left", fontsize=9, fontweight="bold")
axs[0, 1].legend(frameon=False, fontsize=6.5, ncol=4, loc="lower right")
fig.savefig(ROOT / "figures" / "fig_throwing.png"); print("written")
