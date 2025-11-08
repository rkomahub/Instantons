import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------- helpers --------
def read_params(hpp_path="src/parameters.hpp"):
    if not os.path.exists(hpp_path):
        return {}
    txt = open(hpp_path, "r").read()
    def grab(name, cast=float):
        m = re.search(rf"{name}\s*=\s*([0-9.eE+-]+)", txt)
        return cast(m.group(1)) if m else None
    return {
        "N":   grab("N", int),
        "a":   grab("a", float),
        "eta": grab("eta", float),
    }

def load_corr_csv(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path)

    # expected columns from C++ export
    # time,corr_1,dcorr_1,corr_2,d_corr2,corr_3,dcorr_3,
    # log_der_c1,dlog_der_c1,log_der_c2,dlog_der_c2,log_der_c3,dlog_der_c3
    required = [
        "time","corr_1","corr_2","corr_3",
        "log_der_c1","log_der_c2","log_der_c3"
    ]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"{path} missing column '{c}'. Found {list(df.columns)}")
    # errors: some builds use 'd_corr2' instead of 'dcorr_2'
    err_cols_map = {
        "dcorr_1": "dcorr_1",
        "d_corr2": "dcorr_2",
        "dcorr_3": "dcorr_3",
        "dlog_der_c1": "dlog_der_c1",
        "dlog_der_c2": "dlog_der_c2",
        "dlog_der_c3": "dlog_der_c3",
    }
    for old,new in err_cols_map.items():
        if old in df.columns:
            df[new] = df[old]
        elif new not in df.columns:
            df[new] = np.nan  # if no errors provided, keep NaNs

    return df

# -------- read data --------
params = read_params()
eta = params.get("eta", None)

csv_path = "data/correlators_log_derivative.csv"  # uncooled
df = load_corr_csv(csv_path)

tau  = df["time"].values
c1   = df["corr_1"].values   # <x(0)x(τ)>
c2   = df["corr_2"].values   # <x^2(0)x^2(τ)>
c3   = df["corr_3"].values   # <x^3(0)x^3(τ)>
dc1  = df["dcorr_1"].values
dc2  = df["dcorr_2"].values
dc3  = df["dcorr_3"].values

ld1  = df["log_der_c1"].values
ld2  = df["log_der_c2"].values
ld3  = df["log_der_c3"].values
dld1 = df["dlog_der_c1"].values
dld2 = df["dlog_der_c2"].values
dld3 = df["dlog_der_c3"].values

# -------- plot (Fig. 4a and 4b) --------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), sharex=True)

# (a) correlators
ax1.plot(tau, c1, label=r"$\langle x(0)x(\tau)\rangle$", color="C0")
ax1.plot(tau, c2, label=r"$\langle x^2(0)x^2(\tau)\rangle$", color="C1")
ax1.plot(tau, c3, label=r"$\langle x^3(0)x^3(\tau)\rangle$", color="C2")
# optional error bands if present (non-NaN)
if np.any(np.isfinite(dc1)):
    ax1.fill_between(tau, c1-dc1, c1+dc1, color="C0", alpha=0.2, linewidth=0)
if np.any(np.isfinite(dc2)):
    ax1.fill_between(tau, c2-dc2, c2+dc2, color="C1", alpha=0.2, linewidth=0)
if np.any(np.isfinite(dc3)):
    ax1.fill_between(tau, c3-dc3, c3+dc3, color="C2", alpha=0.2, linewidth=0)

ax1.set_xlabel(r"$\tau$")
ax1.set_ylabel(r"$\Pi(\tau)$")
title = r"Correlators $\langle O(0)O(\tau)\rangle$"
if eta is not None:
    title += rf"  ($\eta={eta}$)"
ax1.set_title(title)
ax1.grid(True, alpha=0.3)
ax1.legend()

# (b) log-derivatives
ax2.plot(tau, ld1, label=r"$\frac{d}{d\tau}\log\langle x(0)x(\tau)\rangle$", color="C0")
ax2.plot(tau, ld2, label=r"$\frac{d}{d\tau}\log\langle x^2(0)x^2(\tau)\rangle$", color="C1")
ax2.plot(tau, ld3, label=r"$\frac{d}{d\tau}\log\langle x^3(0)x^3(\tau)\rangle$", color="C2")
# optional error bands
if np.any(np.isfinite(dld1)):
    ax2.fill_between(tau, ld1-dld1, ld1+dld1, color="C0", alpha=0.2, linewidth=0)
if np.any(np.isfinite(dld2)):
    ax2.fill_between(tau, ld2-dld2, ld2+dld2, color="C1", alpha=0.2, linewidth=0)
if np.any(np.isfinite(dld3)):
    ax2.fill_between(tau, ld3-dld3, color="C2", alpha=0.2, linewidth=0, y2=ld3+dld3)

ax2.set_xlabel(r"$\tau$")
ax2.set_ylabel(r"log-derivative")
ax2.set_title(r"Logarithmic derivatives")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()