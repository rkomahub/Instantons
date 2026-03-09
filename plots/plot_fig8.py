import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import eigh

# -----------------
# Inputs
# -----------------
FILE_BLUE = "data/fig8_cooling10.csv"   # MC + 10 cooling sweeps
FILE_RED  = "data/fig8_qmidens.csv"     # switching (non-Gaussian corrected density)

# -----------------
# Physics helpers
# -----------------
def S0(eta):
    return 4.0 * eta**3 / 3.0

def n_1loop(eta):
    return 8.0 * eta**2.5 * np.sqrt(2.0 / np.pi) * np.exp(-S0(eta))

def n_2loop(eta):
    return n_1loop(eta) * np.exp(-(71.0 / 72.0) / S0(eta))

def V_double_well(x, eta):
    return (x**2 - eta**2)**2

# -----------------
# Exact diagonalization in HO basis
# -----------------
def ho_basis_matrix(Nbasis, omega0, eta):
    H = np.zeros((Nbasis, Nbasis))
    c = 1.0 / np.sqrt(omega0)
    A = 1.0
    B = -2.0 * eta**2 - omega0**2 / 4.0
    C = eta**4
    c4 = c**4
    c2 = c**2

    for n in range(Nbasis):
        H[n, n] = 3*A*c4*((n+1)**2 + n**2) + B*c2*(2*n + 1) + omega0*(n + 0.5) + C
        if n+2 < Nbasis:
            H[n, n+2] = A*c4*(4*n+6)*np.sqrt((n+1)*(n+2)) + B*c2*np.sqrt((n+1)*(n+2))
            H[n+2, n] = H[n, n+2]
        if n+4 < Nbasis:
            H[n, n+4] = c4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
            H[n+4, n] = H[n, n+4]
    return H

def deltaE_over_2(eta, Nbasis=60):
    omega0 = 4.0 * eta
    H = ho_basis_matrix(Nbasis, omega0, eta)
    evals, _ = eigh(H)
    evals = np.sort(evals)
    return 0.5 * (evals[1] - evals[0])

# -----------------
# Load data
# -----------------
blue = pd.read_csv(FILE_BLUE)
# expected cols: eta,density_mean,density_err
red = None
try:
    red = pd.read_csv(FILE_RED)  # expected cols: eta,density_ng
except FileNotFoundError:
    red = None

etas_blue = blue["eta"].to_numpy()
dens_blue = blue["density_mean"].to_numpy()
err_blue  = blue["density_err"].to_numpy()

# -----------------
# Build smooth theory curves
# -----------------
eta_grid = np.linspace(0.1, 2.0, 200)

n1 = n_1loop(eta_grid)
n2 = n_2loop(eta_grid)

# ED curve
dE2 = np.array([deltaE_over_2(e) for e in eta_grid])

# -----------------
# Plot
# -----------------
plt.figure(figsize=(7.5, 5.2))

# Blue: MC+cool10
plt.errorbar(etas_blue, dens_blue, yerr=err_blue, fmt="o", ms=4, capsize=3,
             label="MC + 10 cooling sweeps")

# Red: non-Gaussian switching
if red is not None and {"eta","density_ng_mean","density_ng_err"}.issubset(red.columns):
    plt.errorbar(red["eta"], red["density_ng_mean"],
                 yerr=red["density_ng_err"],
                 fmt="s", ms=4, capsize=3,
                 label="Switching (non-Gaussian)")

# Green: 1-loop and 2-loop
plt.plot(eta_grid, n1, "-", label="1-loop semiclassical")
plt.plot(eta_grid, n2, "--", label="2-loop semiclassical")

# Black: ΔE/2
plt.plot(eta_grid, dE2, "-", label=r"$\Delta E(\eta)/2$ (ED)")

plt.xlabel(r"$\eta$")
plt.ylabel(r"instanton density  $n_{I{+}A} = N_{I{+}A}/\beta$")
plt.title("Fig.8: Instanton density vs η")
plt.xlim(0, 2)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()