#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("data")
PLOTS = Path("figures")
PLOTS.mkdir(exist_ok=True)

def main():
    rilm = pd.read_csv(DATA/"fig12_rilm_path.csv")
    heated = pd.read_csv(DATA/"fig12_heated_rilm_path.csv")

    fig, ax = plt.subplots(figsize=(8,5))

    ax.plot(rilm["tau"], rilm["x"],
            label="RILM (classical)",
            color="black", linewidth=2.0)

    ax.plot(heated["tau"], heated["x"],
            label="Heated RILM (10 sweeps)",
            color="mediumseagreen", linewidth=1.2)

    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$x(\tau)$")
    ax.set_title(r"Fig.12: RILM vs Heated RILM")
    ax.legend()
    ax.grid(True)

    fig.tight_layout()
    fig.savefig(PLOTS/"Figure_12.png", dpi=250)

if __name__ == "__main__":
    main()