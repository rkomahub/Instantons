#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from pathlib import Path

# =============================
# Parameters
# =============================
eta = 1.4
Nbasis = 40
omega0 = 4.0 * eta

DATA = Path("data")
FIGS = Path("figures")
FIGS.mkdir(exist_ok=True)

# =============================
# Harmonic oscillator basis Hamiltonian
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
        H[n, n] = 3*A*c4*((n+1)**2 + n**2) + B*c2*(2*n + 1) + omega0*(n + 0.5) + C
        if n+2 < Nbasis:
            H[n, n+2] = A*c4*(4*n+6)*np.sqrt((n+1)*(n+2)) + B*c2*np.sqrt((n+1)*(n+2))
            H[n+2, n] = H[n, n+2]
        if n+4 < Nbasis:
            H[n, n+4] = c4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
            H[n+4, n] = H[n, n+4]
    return H

# =============================
# Operator matrices in HO basis
# =============================
def x_operator(Nbasis, omega0):
    X = np.zeros((Nbasis, Nbasis))
    for n in range(Nbasis-1):
        val = np.sqrt((n+1)/(2*omega0))
        X[n, n+1] = val
        X[n+1, n] = val
    return X

def x_power_operator(Nbasis, omega0, p):
    X = x_operator(Nbasis, omega0)
    Op = X.copy()
    for _ in range(p-1):
        Op = Op @ X
    return Op

# =============================
# Exact correlator
# =============================
def exact_correlator_finite_beta(E, V, Op, tau, beta):
    C = np.zeros_like(tau)
    E0 = E[0]
    for n in range(len(E)):
        amp = V[:,0] @ (Op @ V[:,n])
        w = np.abs(amp)**2
        dE = E[n] - E0
        C += w * (np.exp(-dE * tau) + np.exp(-dE * (beta - tau)))
    return C

def log_derivative(tau, C):
    C = np.maximum(C, 1e-300)
    dt = tau[1] - tau[0]
    out = np.full_like(C, np.nan)
    out[1:-1] = -(np.log(C[2:]) - np.log(C[:-2]))/(2*dt)
    return out

# =============================
# Effective mass from correlator ratio
# =============================
def meff_ratio(tau, C):
    C = np.asarray(C, float)
    tau = np.asarray(tau, float)
    out = np.full_like(C, np.nan)
    dt = tau[1] - tau[0]

    # valid where both points are positive
    ok = (C[:-1] > 0.0) & (C[1:] > 0.0)
    out[:-1][ok] = np.log(C[:-1][ok] / C[1:][ok]) / dt
    return out

def meff_ratio_with_err(tau, C, dC):
    C = np.asarray(C, float)
    dC = np.asarray(dC, float)
    tau = np.asarray(tau, float)

    dt = tau[1] - tau[0]
    m  = np.full_like(C, np.nan)
    dm = np.full_like(C, np.nan)

    ok = (C[:-1] > 0.0) & (C[1:] > 0.0)

    m[:-1][ok] = np.log(C[:-1][ok] / C[1:][ok]) / dt

    rel2 = (dC[:-1][ok] / C[:-1][ok])**2 + (dC[1:][ok] / C[1:][ok])**2
    dm[:-1][ok] = np.sqrt(rel2) / dt
    return m, dm

# =============================
# Load MC data
# =============================
df = pd.read_csv(DATA/"fig4_quantum_correlators.csv")
tau = df["tau"].values

C1_mc = df["C1"].values
C2_mc = df["C2conn"].values
C3_mc = df["C3"].values

dC1 = df["C1_err"].values
dC2 = df["C2conn_err"].values
dC3 = df["C3_err"].values

# =============================
# Exact spectral computation
# =============================
H = ho_basis_matrix(Nbasis, omega0, eta)
E, V = eigh(H)

X1 = x_power_operator(Nbasis, omega0, 1)
X2 = x_power_operator(Nbasis, omega0, 2)
X3 = x_power_operator(Nbasis, omega0, 3)

beta = 800 * 0.05  # or read from params; equals 40.0

C1_exact = exact_correlator_finite_beta(E, V, X1, tau, beta)
C2_exact = exact_correlator_finite_beta(E, V, X2, tau, beta)
C3_exact = exact_correlator_finite_beta(E, V, X3, tau, beta)

# Connected substract mean value for C2
mean_x2 = V[:,0] @ (X2 @ V[:,0])
C2_exact = C2_exact - mean_x2**2

m1_mc, dm1_mc = meff_ratio_with_err(tau, C1_mc, dC1)
m2_mc, dm2_mc = meff_ratio_with_err(tau, C2_mc, dC2)
m3_mc, dm3_mc = meff_ratio_with_err(tau, C3_mc, dC3)

m1_ex = meff_ratio(tau, C1_exact)
m2_ex = meff_ratio(tau, C2_exact)
m3_ex = meff_ratio(tau, C3_exact)

# =============================
# Plot
# =============================
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5), sharex=True)

# Panel (a): correlators (Exact)
ax1.plot(tau, C1_exact, 'k-',  lw=2, label=r"Exact $\langle x(0)x(\tau)\rangle$")
ax1.plot(tau, C2_exact, 'k--', lw=2, label=r"Exact $\langle x^2(0)x^2(\tau)\rangle_{\rm conn}$")
ax1.plot(tau, C3_exact, 'k:',  lw=2, label=r"Exact $\langle x^3(0)x^3(\tau)\rangle$")

ax1.errorbar(tau, C1_mc, yerr=dC1, fmt='o', ms=3, label=r"MC $\langle x x\rangle$")
ax1.errorbar(tau, C2_mc, yerr=dC2, fmt='o', ms=3, label=r"MC $\langle x^2 x^2\rangle_{\rm conn}$")
ax1.errorbar(tau, C3_mc, yerr=dC3, fmt='o', ms=3, label=r"MC $\langle x^3 x^3\rangle$")

ax1.set_xlim(0.0, 1.5)
ax1.set_ylim(0.0, 8.0)
ax1.set_xlabel(r"$\tau$")
ax1.set_ylabel(r"$C(\tau)$")
ax1.grid(True, alpha=0.3)
ax1.legend()

# Panel (b): effective masses (Exact)
ax2.plot(tau, m1_ex, 'k-',  lw=2, label=r"Exact $m_{\rm eff}^{(x)}$")
ax2.plot(tau, m2_ex, 'k--', lw=2, label=r"Exact $m_{\rm eff}^{(x^2)}$")
ax2.plot(tau, m3_ex, 'k:',  lw=2, label=r"Exact $m_{\rm eff}^{(x^3)}$")

ax2.errorbar(tau, m1_mc, yerr=dm1_mc, fmt='o', ms=3, label=r"MC $x$")
ax2.errorbar(tau, m2_mc, yerr=dm2_mc, fmt='o', ms=3, label=r"MC $x^2$")
ax2.errorbar(tau, m3_mc, yerr=dm3_mc, fmt='o', ms=3, label=r"MC $x^3$")

ax2.set_xlim(0.0, 1.5)
ax2.set_ylim(0.0, 8.0)   # this is what Schäfer effectively shows
ax2.set_xlabel(r"$\tau$")
ax2.set_ylabel(r"$m_{\mathrm{eff}}(\tau)$")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig(FIGS/"Figure_4.png", dpi=250)
plt.show()