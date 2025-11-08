import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------- helpers ---------
def read_params(hpp_path="src/parameters.hpp"):
    if not os.path.exists(hpp_path):
        return None
    txt = open(hpp_path, "r").read()
    def grab(name, cast=float):
        m = re.search(rf"{name}\s*=\s*([0-9.eE+-]+)", txt)
        return cast(m.group(1)) if m else None
    return {
        "N":   grab("N", int),
        "a":   grab("a", float),
        "eta": grab("eta", float),
    }

def load_csv_required(path, cols):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")
    df = pd.read_csv(path)
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"{path} missing column '{c}'. Found: {list(df.columns)}")
    return df

# --------- load C++ data ---------
eig_path  = "data/eigenvalues_eigenfunctions.csv"      # from your diagonalizer main
hist_path = "data/probability_histogram.csv"           # from your MC Metropolis main

eig = load_csv_required(eig_path,
                        ["positions","eigenvalues","eigenfunction0","eigenfunction1","eigenfunction2","eigenfunction3"])
hist = load_csv_required(hist_path, ["probability"])

params = read_params()  # optional; used for a nicer title
eta = params["eta"] if params and params["eta"] is not None else "?"

# exact |psi0|^2 from C++ output
x_grid = eig["positions"].values
psi0   = eig["eigenfunction0"].values
psi0_sq = np.abs(psi0)**2

# re-normalize |psi0|^2 just in case (should already be normalized, but this is safe)
norm = np.trapz(psi0_sq, x_grid)
if norm > 0:
    psi0_sq = psi0_sq / norm

# MC histogram counts → probability density
counts = hist["probability"].values.astype(float)

# Use the eigenfunction x-range to define the histogram bin edges/centers,
# assuming the MC histogram was filled over the same interval.
xmin, xmax = np.min(x_grid), np.max(x_grid)
nbins = len(counts)
bin_edges = np.linspace(xmin, xmax, nbins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
bin_width = bin_edges[1] - bin_edges[0]

# If counts are raw, convert to PDF. If they’re already normalized, this won’t break things.
total_counts = counts.sum()
if total_counts > 0 and bin_width > 0:
    pdf = counts / (total_counts * bin_width)
else:
    pdf = counts  # fallback

# --------- plot ---------
plt.figure(figsize=(7.2,5.2))

# MC probability density (step)
plt.step(bin_centers, pdf, where="mid", color="tab:blue", label="Monte Carlo", alpha=0.9)

# Exact |psi0|^2 from C++ eigenfile (line)
plt.plot(x_grid, psi0_sq, color="black", linewidth=2.0, label=r"exact $|\psi_0(x)|^2$")

plt.xlabel(r"$x$")
plt.ylabel("Probability density")
title = r"Probability distribution in double well"
if eta != "?":
    title += rf"  ($\eta={eta}$)"
plt.title(title)
plt.xlim(xmin, xmax)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()