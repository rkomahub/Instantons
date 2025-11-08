#!/usr/bin/env python3
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

PARAMS = Path("src/parameters.hpp")
DATA   = Path("data")

def read_params(hpp=PARAMS):
    txt = Path(hpp).read_text()
    def grab(name, cast=float):
        m = re.search(rf"{name}\s*=\s*([0-9.eE+-]+)", txt)
        if not m:
            raise ValueError(f"Parameter '{name}' not found in {hpp}")
        return cast(m.group(1))
    N   = grab("N", int)
    a   = grab("a", float)
    eta = grab("eta", float)
    return N, a, eta, N*a

def load_csv_or_die(path, expected_cols=None):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}\nRun your C++ binary first to generate it.")
    df = pd.read_csv(path)
    if expected_cols and not all(c in df.columns for c in expected_cols):
        raise ValueError(f"{path} must contain columns {expected_cols}, got {list(df.columns)}")
    return df

def log_derivative(tau, C, subtract_const=False):
    """
    Discrete log-derivative like Schäfer Fig. 4b/6:
      d/dτ log C(τ) ≈ (log C(τ+Δτ) - log C(τ)) / Δτ
    If subtract_const=True (even powers), subtract C(τ_max) before logging.
    """
    Cuse = C.copy()
    if subtract_const:
        Cuse = Cuse - Cuse.iloc[-1]
    # guard against non-positive values
    eps = 1e-14
    Cuse = Cuse.clip(lower=eps)
    d = (Cuse.shift(-1).apply(math.log) - Cuse.apply(math.log)) / (tau.shift(-1) - tau)
    return d.iloc[:-1], tau.iloc[:-1]

def main():
    N, a, eta, beta = read_params()

    # --- Fig. 2: paths ---
    q_path = load_csv_or_die(DATA/"quantum_path.csv", ["tau","x"])
    c_path = load_csv_or_die(DATA/"cooled_path.csv",  ["tau","x"])

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(q_path["tau"], q_path["x"], label="Monte Carlo", color="black")
    ax.plot(c_path["tau"], c_path["x"], label="cooled",      color="green")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$x$")
    ax.set_title(rf"Typical Euclidean path  (η={eta}, N={N}, a={a}, β={beta:.2f})")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(DATA/"fig2_paths.png", dpi=200)

    # --- Fig. 6-style: correlators and log-derivatives ---
    q_corr = load_csv_or_die(DATA/"quantum_correlator.csv", ["tau","C(tau)"])
    c_corr = load_csv_or_die(DATA/"cooled_correlator.csv",  ["tau","C(tau)"])

    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
    # (a) correlators
    ax1.plot(q_corr["tau"], q_corr["C(tau)"], label="⟨x(0)x(τ)⟩ (MC)", color="black")
    ax1.plot(c_corr["tau"], c_corr["C(tau)"], label="⟨x(0)x(τ)⟩ (cooled)", color="green")
    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"$\Pi(\tau)$")
    ax1.set_title("Two-point correlator")
    ax1.legend()
    ax1.grid(True)

    # (b) log-derivative ~ energy gap extractor
    dlog_q, tau_q = log_derivative(q_corr["tau"], q_corr["C(tau)"], subtract_const=False)
    dlog_c, tau_c = log_derivative(c_corr["tau"], c_corr["C(tau)"], subtract_const=False)
    ax2.plot(tau_q, dlog_q, label="MC", color="black")
    ax2.plot(tau_c, dlog_c, label="cooled", color="green")
    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$\frac{d}{d\tau}\log \Pi(\tau)$")
    ax2.set_title("Log-derivative (gap → plateau)")
    ax2.legend()
    ax2.grid(True)

    fig2.tight_layout()
    fig2.savefig(DATA/"fig6_correlators_logder.png", dpi=200)

    print("Saved:")
    print(" - data/fig2_paths.png")
    print(" - data/fig6_correlators_logder.png")

if __name__ == "__main__":
    main()
