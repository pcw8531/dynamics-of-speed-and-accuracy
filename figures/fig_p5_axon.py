"""Mean score by play in the Axon data with the fit of Equation 18 (analysis/p5_mean_curve.csv)."""
import logging
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
ROOT = Path(__file__).resolve().parent.parent
c = pd.read_csv(ROOT / "analysis" / "p5_mean_curve.csv")
plt.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "axes.linewidth": 0.7})
fig, ax = plt.subplots(figsize=(3.3, 2.4)); fig.subplots_adjust(left=0.2, bottom=0.18, right=0.97, top=0.95)
ax.plot(c.play, c.mean_score / 1000, "o", ms=3, color="black", mec="white", mew=0.4, label="mean score")
ax.plot(c.play, c.fit_eq18 / 1000, color="black", lw=1.0, label="Equation 18")
ax.set_xlabel("play"); ax.set_ylabel("score (thousands)"); ax.legend(frameon=False, loc="lower right", fontsize=7.5)
fig.savefig(ROOT / "figures" / "fig_p5_axon.png"); print("written")
