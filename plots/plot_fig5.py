#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from pathlib import Path

# =============================
# Parameters (must match C++)
# =============================
eta = 1.4
Nbasis = 120
omega0 = 4.0 * eta

DATA = Path("data")
FIGS = Path("figures")
FIGS.mkdir(exist_ok=True)

# =============================
# Hamiltonian (exact spectrum)
# =============================
def ho_basis_matrix(Nbasis, omega0, eta):
    H = np.zeros((Nbasis, Nbasis))
    c = 1.0 / np.sqrt(omega0)
    A = 1.0
    B = -2.0 * eta**2 - omega0**2 / 4.0
    C = eta**4
    c4 = c**4
    c2 = c**2

    for n in range(Nbasis):
        H[n, n] = (
            3.0 * A * c4 * ((n + 1)**2 + n**2)
            + B * c2 * (2 * n + 1)
            + omega0 * (n + 0.5)
            + C
        )

        if n + 2 < Nbasis:
            H[n, n + 2] = (
                A * c4 * (4 * n + 6) * np.sqrt((n + 1) * (n + 2))
                + B * c2 * np.sqrt((n + 1) * (n + 2))
            )
            H[n + 2, n] = H[n, n + 2]

        if n + 4 < Nbasis:
            H[n, n + 4] = c4 * np.sqrt((n + 1) * (n + 2) * (n + 3) * (n + 4))
            H[n + 4, n] = H[n, n + 4]

    return H

H = ho_basis_matrix(Nbasis, omega0, eta)
E, _ = eigh(H)

# =============================
# Exact free energy
# =============================
def free_energy_exact(beta):
    Z = np.sum(np.exp(-beta * E))
    return -(1.0 / beta) * np.log(Z)

# =============================
# Load Monte Carlo data
# =============================
df = pd.read_csv(DATA / "fig5_free_energy.csv")
df = df.sort_values("T")

beta_mc = df["beta"].to_numpy()
T_mc = df["T"].to_numpy()

# C++ stores the true thermodynamic free energy F
F_mc_mean = df["F_mean"].to_numpy()
F_mc_err = df["F_err"].to_numpy()

# To match Schäfer's plotted convention, plot -F
F_mc_plot = F_mc_mean

# =============================
# Exact curve
# =============================
beta_exact = np.linspace(1, beta_mc.max(), 400)
T_exact = 1.0 / beta_exact
F_exact = np.array([free_energy_exact(b) for b in beta_exact])

# To match Schäfer's plotted convention, plot -F
F_exact_plot = F_exact

order = np.argsort(T_exact)
T_exact = T_exact[order]
F_exact_plot = F_exact_plot[order]

# =============================
# Plot
# =============================
plt.figure(figsize=(6, 5))

plt.plot(T_exact, F_exact_plot, lw=2, label="Exact", color="orange", linestyle="--")

plt.errorbar(
    T_mc,
    F_mc_plot,
    yerr=F_mc_err,
    fmt="o",
    capsize=3,
    label="MC (switching)",
    color="orange"
)

plt.xlabel(r"$T$")
plt.ylabel(r"$F(T)$")
plt.xscale("log")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGS / "Figure_5.png", dpi=300)
plt.show()