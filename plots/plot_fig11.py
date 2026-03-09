# plots/plot_fig11.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ----------------------------
# Inputs
# ----------------------------
eta = 1.4
csv_path = "data/fig11_gaussian_effective_potential.csv"
out_path = "plots/fig11_gaussian_effective_potential_3D.png"

df = pd.read_csv(csv_path)
tau = df["tau"].to_numpy()
xI_tau = df["xI"].to_numpy()

# ----------------------------
# Build smooth surface V_gauss(x, tau)
# ----------------------------
# If you want the surface smoother, increase Nx, Nt
Nt = 260
Nx = 260

tau_s = np.linspace(tau.min(), tau.max(), Nt)
x_s = np.linspace(-2.1 * eta, 2.1 * eta, Nx)

T, X = np.meshgrid(tau_s, x_s)

# Smooth instanton line via interpolation
xI = np.interp(tau_s, tau, xI_tau)              # shape (Nt,)
xI_mesh = np.broadcast_to(xI, T.shape)          # shape (Nx, Nt)

# V''(x) for V(x) = (x^2 - eta^2)^2
Vpp = 12.0 * xI_mesh**2 - 4.0 * eta**2

# Gaussian effective potential around the valley x = xI(tau)
V = 0.5 * Vpp * (X - xI_mesh)**2

# ----------------------------
# Figure styling
# ----------------------------
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "figure.dpi": 140,
})

fig = plt.figure(figsize=(9.2, 6.2))
ax = fig.add_subplot(111, projection="3d")

# Smooth surface: no edge lines, antialiased, good colormap
surf = ax.plot_surface(
    T, X, V,
    rcount=200, ccount=200,
    linewidth=0,
    antialiased=True,
    shade=True,
    cmap="viridis"
)

# Colorbar
cbar = fig.colorbar(surf, ax=ax, shrink=0.72, pad=0.08, aspect=18)
cbar.set_label(r"$V_{\mathrm{gauss}}(x,\tau)$")

# ----------------------------
# Overlay the instanton trajectory and its projection
# ----------------------------
# Valley line (x = xI(tau), V=0)
tau_line = tau_s
x_line = xI
V_line = np.zeros_like(tau_line)

ax.plot(tau_line, x_line, V_line, linewidth=2.2)

# Project trajectory onto the base plane z = zmin
zmin = 0.0
ax.plot(tau_line, x_line, zmin * np.ones_like(tau_line), linewidth=1.6)

# Optionally add a faint reference line x=0 on the base plane
ax.plot(tau_line, np.zeros_like(tau_line), zmin * np.ones_like(tau_line), linewidth=1.0)

# ----------------------------
# Axes labels and limits
# ----------------------------
ax.set_xlabel(r"$\tau$")
ax.set_ylabel(r"$x$")
ax.set_zlabel(r"$V_{\mathrm{gauss}}(x,\tau)$")

ax.set_xlim(tau_s.min(), tau_s.max())
ax.set_ylim(x_s.min(), x_s.max())
ax.set_zlim(zmin, np.percentile(V, 99.5))

# Nice ticks
ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
ax.zaxis.set_major_locator(ticker.MaxNLocator(6))

# View angle tuned for “article-like” look
ax.view_init(elev=25, azim=-55)

# Clean layout
fig.tight_layout()
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.show()