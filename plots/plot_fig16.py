import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

eta = 1.4

# -----------------------------
# Load separation data
# -----------------------------
df = pd.read_csv("data/fig16_separations.csv")

mc = df[df["source"] == "mc"].copy()
ref = df[df["source"] == "ref"].copy()

# restrict to displayed range
mc = mc[mc["sep"] <= 4.0]
ref = ref[ref["sep"] <= 4.0]

print("len(mc)  =", len(mc))
print("len(ref) =", len(ref))
print(df.head())

# -----------------------------
# Histogram settings
# -----------------------------
bins = np.linspace(0.0, 4.0, 41)
centers = 0.5 * (bins[:-1] + bins[1:])

hist_mc, _ = np.histogram(mc["sep"], bins=bins, density=False)
hist_ref, _ = np.histogram(ref["sep"], bins=bins, density=False)

# -----------------------------
# Streamline IA interaction theory
# -----------------------------
ia = pd.read_csv("data/fig14_streamline_tauZ.csv")

tau_th = ia["tauZ"].to_numpy()
Sint_over_S0 = ia["Sint_over_S0"].to_numpy()

# keep only the visible range
mask = (tau_th >= 0.0) & (tau_th <= 4.0)
tau_th = tau_th[mask]
Sint_over_S0 = Sint_over_S0[mask]

# sort by tau just to be safe
order = np.argsort(tau_th)
tau_th = tau_th[order]
Sint_over_S0 = Sint_over_S0[order]

# classical instanton action
S0 = 4.0 * eta**3 / 3.0

# interaction action
Sint = Sint_over_S0 * S0

# Boltzmann weight
W = np.exp(-Sint)

# normalize as probability density
W /= np.trapezoid(W, tau_th)

# convert to histogram count scale
bin_width = bins[1] - bins[0]
Ntot = len(mc["sep"])
n_theory = W * Ntot * bin_width

# -----------------------------
# Diagnostics
# -----------------------------
print("mc sep min/max =", mc["sep"].min(), mc["sep"].max())
print("ref sep min/max =", ref["sep"].min(), ref["sep"].max())
print("hist_mc max =", hist_mc.max(), "sum =", hist_mc.sum())
print("hist_ref max =", hist_ref.max(), "sum =", hist_ref.sum())
print("n_theory max =", n_theory.max())

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(7, 5))

plt.step(centers, hist_mc, where="mid", color="black", label="MC (cooled)")
plt.step(centers, hist_ref, where="mid", color="blue", label="Random reference")
plt.plot(tau_th, n_theory, "r-", lw=2, label="IA interaction theory (streamline)")

plt.xlim(0, 4)
plt.ylim(0, 1.1 * max(hist_mc.max(), hist_ref.max(), n_theory.max()))
plt.xlabel(r"$\Delta \tau$")
plt.ylabel(r"$n_{IA}(\Delta \tau)$")
plt.title("Instanton–Anti-Instanton separation distribution")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()