#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("data")
PLOTS = Path("figures")
PLOTS.mkdir(exist_ok=True)

def main():
    q_path = pd.read_csv(DATA/"quantum_path.csv")
    c_path = pd.read_csv(DATA/"cooled_path.csv")

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(q_path["tau"], q_path["x"], label="Monte Carlo", color="black", linewidth=1.0)
    ax.plot(c_path["tau"], c_path["x"], label="Cooled (200 sweeps)", color="mediumseagreen", linewidth=2.0)

    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$x(\tau)$")
    ax.set_title(r"Fig. 2: Typical Euclidean path ($\eta=1.4,\ a=0.05,\ N_\tau=800$)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(PLOTS/"Figure_2.png", dpi=250)

if __name__ == "__main__":
    main()