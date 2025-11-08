import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigh

# -----------------
# Parameters
# -----------------
eta = 1.4
x_min, x_max = -2.2, 2.2
Nx = 400  # grid points for plotting potential

# -----------------
# Potential function
# -----------------
def V_double_well(x, eta):
    return (x**2 - eta**2)**2

# -----------------
# Build Hamiltonian in harmonic oscillator basis
# -----------------
def ho_basis_matrix(Nbasis, omega0, eta):
    """
    Construct Hamiltonian matrix H in harmonic oscillator basis
    following Schäfer's eqs. (5)-(7).
    """
    # Matrix
    H = np.zeros((Nbasis, Nbasis))
    c = 1.0 / np.sqrt(omega0)
    A = 1.0
    B = -2.0 * eta**2 - omega0**2 / 4.0
    C = eta**4
    c4 = c**4
    c2 = c**2

    for n in range(Nbasis):
        # diagonal
        H[n, n] = 3*A*c4*((n+1)**2 + n**2) + B*c2*(2*n + 1) + omega0*(n + 0.5) + C
        # n → n+2
        if n+2 < Nbasis:
            H[n, n+2] = A*c4*(4*n+6)*np.sqrt((n+1)*(n+2)) + B*c2*np.sqrt((n+1)*(n+2))
            H[n+2, n] = H[n, n+2]
        # n → n+4
        if n+4 < Nbasis:
            H[n, n+4] = c4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
            H[n+4, n] = H[n, n+4]

    return H

# -----------------
# Figure 1a — Potential with energy levels
# -----------------
x_vals = np.linspace(x_min, x_max, Nx)
V_vals = V_double_well(x_vals, eta)

# Diagonalize for fixed eta
Nbasis = 40
omega0 = 4.0 * eta
Hmat = ho_basis_matrix(Nbasis, omega0, eta)
evals, _ = eigh(Hmat)
evals = np.sort(evals)[:4]  # first four energies

plt.figure(figsize=(6, 5))
plt.plot(x_vals, V_vals, 'k-', label=r'$V(x)$')
colors = ['r', 'g', 'b', 'm']
for i, E in enumerate(evals):
    plt.hlines(E, x_min, x_max, colors=colors[i], linestyles='--', label=f'E{i}')
plt.xlabel(r'$x$')
plt.ylabel(r'$V(x)$')
plt.ylim(0, max(V_vals)*1.1)
plt.title(f'Double well potential, η={eta}')
plt.legend()
plt.grid(True)

# -----------------
# Figure 1b — Spectrum vs η
# -----------------
eta_values = np.linspace(0.5, 2.0, 30)
Nbasis = 40
n_states = 6
energies_vs_eta = np.zeros((len(eta_values), n_states))

for idx, eta_val in enumerate(eta_values):
    omega0 = 4.0 * eta_val
    Hmat = ho_basis_matrix(Nbasis, omega0, eta_val)
    evals, _ = eigh(Hmat)
    evals = np.sort(evals)[:n_states]
    energies_vs_eta[idx, :] = evals

plt.figure(figsize=(6, 5))
for n in range(n_states):
    plt.plot(eta_values, energies_vs_eta[:, n], label=f'E{n}')
plt.xlabel(r'$\eta$')
plt.ylabel('Energy')
plt.title('Spectrum vs η')
plt.legend()
plt.grid(True)

plt.show()