#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from numpy.linalg import eigh

DATA = Path("data")
FIGS = Path("figures")
FIGS.mkdir(exist_ok=True)

# --- parameters (must match C++) ---
N = 800
a = 0.05
eta = 1.4

# -------------------------------------
# Load ensemble position samples
# -------------------------------------
df = pd.read_csv(DATA / "ensemble_quantum_positions.csv")
x_samples = df["x"].values

# Histogram (Monte Carlo)
nbins = 80
counts, edges = np.histogram(x_samples, bins=nbins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])

# -------------------------------------
# Exact ground state |psi0|^2
# -------------------------------------
# Harmonic oscillator basis truncation
x_grid = np.linspace(-3*eta, 3*eta, 400)
dx = x_grid[1] - x_grid[0]

# Hamiltonian in coordinate space (finite difference)
diag = 2.0 / dx**2 + (x_grid**2 - eta**2)**2
off  = -1.0 / dx**2 * np.ones(len(x_grid)-1)

H = np.diag(diag) + np.diag(off,1) + np.diag(off,-1)

eigvals, eigvecs = eigh(H)
psi0 = eigvecs[:,0]
psi0_sq = np.abs(psi0)**2
psi0_sq /= np.trapezoid(psi0_sq, x_grid)

# -------------------------------------
# Plot
# -------------------------------------
plt.figure(figsize=(7,5))

plt.step(centers, counts, where="mid",
         color="mediumseagreen", label="Monte Carlo")

plt.plot(x_grid, psi0_sq,
         color="black", linewidth=2,
         label=r"Exact $|\psi_0(x)|^2$")

plt.xlabel(r"$x$")
plt.ylabel("Probability density")
plt.title(r"Fig. 3: Probability distribution ($\eta=1.4$)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(FIGS / "Figure_3.png", dpi=250)
plt.savefig(FIGS / "Figure_3.pdf")
plt.show()