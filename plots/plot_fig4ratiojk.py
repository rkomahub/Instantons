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
Nbasis = 80
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
        val = np.sqrt((n+1)/(omega0)) # m=1/2 cancelled out, so just 1/sqrt(omega0) instead of 1/sqrt(2*omega0)
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
# Exact correlator (ground-state dominated, symmetric in tau)
# ============================= 
def exact_correlator_ground_state(E, V, Op, tau, nstates=None):
    """
    Ground-state correlator:
        C(tau) = sum_n |<0|Op|n>|^2 * exp(-(E_n - E_0) * tau)

    Parameters
    ----------
    E : (Nbasis,) eigenvalues
    V : (Nbasis, Nbasis) eigenvectors, columns V[:,n] are |n>
    Op : (Nbasis, Nbasis) operator matrix
    tau : array of Euclidean times
    nstates : optional int, number of states to include in the sum

    Returns
    -------
    C : array same shape as tau
    """
    tau = np.asarray(tau, float)
    C = np.zeros_like(tau)

    E0 = E[0]
    nmax = len(E) if nstates is None else min(nstates, len(E))

    v0 = V[:, 0]
    Op_v = Op @ V[:, :nmax]               # (Nbasis, nmax)
    amps = v0 @ Op_v                      # (nmax,)
    w = np.abs(amps)**2                   # |<0|Op|n>|^2
    dE = E[:nmax] - E0

    # include n=0 term automatically (gives constant = |<0|Op|0>|^2)
    C = (w[:, None] * np.exp(-dE[:, None] * tau[None, :])).sum(axis=0)
    return C

# =============================
# Effective mass from correlator ratio
# =============================

def meff_ratio(tau, C):
    C = np.asarray(C, float)
    tau = np.asarray(tau, float)
    out = np.full_like(C, np.nan)
    dt = tau[1] - tau[0]

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

def meff_ratio_from_C(Cmean, dt):
    """
    Cmean: (N,) mean correlator
    returns m: (N,) with last point NaN
    """
    Cmean = np.asarray(Cmean, float)
    m = np.full_like(Cmean, np.nan)
    ok = (Cmean[:-1] > 0.0) & (Cmean[1:] > 0.0)
    m[:-1][ok] = np.log(Cmean[:-1][ok] / Cmean[1:][ok]) / dt
    return m

# =============================
# Jackknife effective mass from correlator ratio
# =============================

def jackknife_meff_ratio(C_trials, tau):
    """
    Jackknife on the nonlinear estimator m_eff = log(C_i/C_{i+1}) / dt,
    using leave-one-out means of C.
    """
    C_trials = np.asarray(C_trials, float)
    tau = np.asarray(tau, float)

    T, N = C_trials.shape
    dt = tau[1] - tau[0]

    sumC = C_trials.sum(axis=0)                 # (N,)
    # leave-one-out means: C_loo[k] = (sumC - C_k)/(T-1)
    C_loo = (sumC[None, :] - C_trials) / (T - 1)  # (T, N)

    m_loo = np.array([meff_ratio_from_C(C_loo[k], dt) for k in range(T)])  # (T, N)

    m_mean = np.nanmean(m_loo, axis=0)
    var = (T - 1) / T * np.nansum((m_loo - m_mean[None, :])**2, axis=0)
    m_err = np.sqrt(var)

    return m_mean, m_err

# =============================
# Load MC data (now includes C2raw and C2conn)
# =============================
df = pd.read_csv(DATA / "fig4_quantum_correlators.csv")
tau = df["tau"].values

C1_mc = df["C1"].values
dC1 = df["C1_err"].values

C2raw_mc = df["C2raw"].values
dC2raw = df["C2raw_err"].values

C2conn_mc = df["C2conn"].values
dC2conn = df["C2conn_err"].values

C3_mc = df["C3"].values
dC3 = df["C3_err"].values

# =============================

dfT = pd.read_csv(DATA / "fig4_trials_correlators.csv")

# Ensure same tau ordering as main file
taus = np.sort(dfT["tau"].unique())
if len(taus) != len(tau) or np.max(np.abs(taus - tau)) > 1e-12:
    raise RuntimeError("tau grid in fig4_trials_correlators.csv does not match fig4_quantum_correlators.csv")

T = int(dfT["trial"].nunique())
N = len(tau)

# Pivot to (T, N) arrays
C1_trials = dfT.pivot(index="trial", columns="tau", values="C1").to_numpy()
C2conn_trials = dfT.pivot(index="trial", columns="tau", values="C2conn").to_numpy()
C3_trials = dfT.pivot(index="trial", columns="tau", values="C3").to_numpy()

# Safety
if C1_trials.shape != (T, N):
    raise RuntimeError(f"Unexpected C1_trials shape {C1_trials.shape}, expected {(T,N)}")

# =============================
# Exact spectral computation
# =============================
H = ho_basis_matrix(Nbasis, omega0, eta)
E, V = eigh(H)

X1 = x_power_operator(Nbasis, omega0, 1)
X2 = x_power_operator(Nbasis, omega0, 2)
X3 = x_power_operator(Nbasis, omega0, 3)

# Use ground-state spectral correlator (T -> 0 form)
C1_exact = exact_correlator_ground_state(E, V, X1, tau, nstates=20)

# Raw and connected x^2 correlators separated
C2_exact_raw = exact_correlator_ground_state(E, V, X2, tau, nstates=30)
mean_x2 = V[:, 0] @ (X2 @ V[:, 0])
C2_exact_conn = C2_exact_raw - mean_x2**2

C3_exact = exact_correlator_ground_state(E, V, X3, tau, nstates=40)  # or None

# =============================
# Jackknife effective masses (ratio)
#   - x: raw correlator
#   - x^2: connected correlator
#   - x^3: raw correlator
# =============================
m1_mc, dm1_mc = jackknife_meff_ratio(C1_trials, tau)
m2_mc, dm2_mc = jackknife_meff_ratio(C2conn_trials, tau) # connected x^2 correlator
m3_mc, dm3_mc = jackknife_meff_ratio(C3_trials, tau)

m1_ex = meff_ratio(tau, C1_exact)
m2_ex = meff_ratio(tau, C2_exact_conn)
m3_ex = meff_ratio(tau, C3_exact)

# =============================
# Gap comparison diagnostic
# =============================

# Exact spectral gaps
Delta1 = E[1] - E[0]
Delta2 = E[2] - E[0]

print("\n=== Exact spectral gaps ===")
print(f"E1 - E0 = {Delta1:.6f}")
print(f"E2 - E0 = {Delta2:.6f}")

# Choose plateau window
tau_min = 1.0
tau_max = 1.5

mask = (tau >= tau_min) & (tau <= tau_max)

def plateau_average(m):
    return np.nanmean(m[mask])

print("\n=== MC plateau estimates ===")
print(f"x      plateau ≈ {plateau_average(m1_mc):.6f}")
print(f"x^3    plateau ≈ {plateau_average(m3_mc):.6f}")
print(f"x^2(conn) plateau ≈ {plateau_average(m2_mc):.6f}")

print("\n=== Exact curve plateau estimates ===")
print(f"x      plateau ≈ {plateau_average(m1_ex):.6f}")
print(f"x^3    plateau ≈ {plateau_average(m3_ex):.6f}")
print(f"x^2(conn) plateau ≈ {plateau_average(m2_ex):.6f}")

# Exact <x^2>
x2_exact = V[:,0] @ (X2 @ V[:,0])

# Estimate <x^2>_MC from large-tau limit of raw correlator
tau_min = 1.0
tau_max = 1.5
mask = (tau >= tau_min) & (tau <= tau_max)

C2_plateau = np.nanmean(C2raw_mc[mask])
x2_mc_est = np.sqrt(C2_plateau)

print("\n=== Exact and Montecarlo estimate <x^2> ===")
print(f"<x^2>_exact = {x2_exact:.6f}")
print(f"<x^2>_MC (from plateau) ≈ {x2_mc_est:.6f}")

# =============================
# Plot
# =============================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

# Panel (a): correlators (raw for x^2)
ax1.plot(tau, C1_exact, 'k-',  lw=2, label=r"Exact $\langle x(0)x(\tau)\rangle$")
ax1.plot(tau, C2_exact_raw, 'k--', lw=2, label=r"Exact $\langle x^2(0)x^2(\tau)\rangle$")
ax1.plot(tau, C3_exact, 'k:',  lw=2, label=r"Exact $\langle x^3(0)x^3(\tau)\rangle$")

ax1.errorbar(tau, C1_mc, yerr=dC1, fmt='o', ms=3, label=r"MC $\langle x x\rangle$")
ax1.errorbar(tau, C2raw_mc, yerr=dC2raw, fmt='o', ms=3, label=r"MC $\langle x^2 x^2\rangle$")
ax1.errorbar(tau, C3_mc, yerr=dC3, fmt='o', ms=3, label=r"MC $\langle x^3 x^3\rangle$")

ax1.set_xlim(0.0, 1.5)
ax1.set_ylim(0.0, 8.0)
ax1.set_xlabel(r"$\tau$")
ax1.set_ylabel(r"$C(\tau)$")
ax1.grid(True, alpha=0.3)
ax1.legend()

# Panel (b): effective masses (ratio, connected for x^2)
ax2.plot(tau, m1_ex, 'k-',  lw=2, label=r"Exact $m_{\rm eff}^{(x)}$")
ax2.plot(tau, m2_ex, 'k--', lw=2, label=r"Exact $m_{\rm eff}^{(x^2)}$ (conn)")
ax2.plot(tau, m3_ex, 'k:',  lw=2, label=r"Exact $m_{\rm eff}^{(x^3)}$")

ax2.errorbar(tau, m1_mc, yerr=dm1_mc, fmt='o', ms=3, label=r"MC $x$")
ax2.errorbar(tau, m2_mc, yerr=dm2_mc, fmt='o', ms=3, label=r"MC $x^2$ (conn)")
ax2.errorbar(tau, m3_mc, yerr=dm3_mc, fmt='o', ms=3, label=r"MC $x^3$")

ax2.axhline(Delta1, color='gray', linestyle='--', alpha=0.6, label=r"$E_1 - E_0$")
ax2.axhline(Delta2, color='gray', linestyle=':',  alpha=0.6, label=r"$E_2 - E_0$")

ax2.set_xlim(0.0, 1.5)
ax2.set_ylim(0.0, 8.0)
ax2.set_xlabel(r"$\tau$")
ax2.set_ylabel(r"$m_{\mathrm{eff}}(\tau)$")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig(FIGS / "Figure_4ratiojk.png", dpi=250)
plt.show()