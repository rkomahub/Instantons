import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import defaultdict

CSV_PATH = "data/fig9_paths.csv"
OUT_PATH = "figures/Figure_9.png"

ALPHA_TO_PLOT = 1.0
SAMPLE_ID_TO_PLOT = 0

data = defaultdict(lambda: {"tau": [], "x": []})

with open(CSV_PATH, "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        sector = int(row["sector"])
        alpha = float(row["alpha"])
        sample_id = int(row["sample_id"])
        tau = float(row["tau"])
        x = float(row["x"])
        data[(sector, alpha, sample_id)]["tau"].append(tau)
        data[(sector, alpha, sample_id)]["x"].append(x)

for key in data:
    tau = np.array(data[key]["tau"])
    x = np.array(data[key]["x"])
    idx = np.argsort(tau)
    data[key] = (tau[idx], x[idx])

fig, ax = plt.subplots(figsize=(7, 5))

sector_color = {0: "tab:blue", 1: "tab:green"}
sector_label = {0: "sector 0", 1: "sector 1"}

for sector in [0, 1]:
    key = (sector, ALPHA_TO_PLOT, SAMPLE_ID_TO_PLOT)
    if key not in data:
        continue

    tau, x = data[key]
    ax.plot(tau, x, color=sector_color[sector], linewidth=1.5,
            label=sector_label[sector])

ax.set_title(rf"Fig.9: Switching paths at $\alpha = {ALPHA_TO_PLOT:g}$")
ax.set_xlabel(r"$\tau$")
ax.set_xlim(0, 35)
ax.set_ylabel(r"$x(\tau)$")
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
print(f"[✓] Saved {OUT_PATH}")
plt.show()