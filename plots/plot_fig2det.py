#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("data")
FIGS = Path("figures")
FIGS.mkdir(exist_ok=True)

# Must match your C++ parameters for the run
A = 0.05
ETA = 1.4


def V(x: np.ndarray) -> np.ndarray:
    """Double well potential V(x) = (x^2 - eta^2)^2."""
    return (x * x - ETA * ETA) ** 2


def dVdx(x: np.ndarray) -> np.ndarray:
    """Derivative V'(x) = 4 x (x^2 - eta^2)."""
    return 4.0 * x * (x * x - ETA * ETA)


def action(x: np.ndarray) -> float:
    """
    Discretized Euclidean action with periodic BC:
      S = sum_i [ (x_i - x_{i-1})^2 / (4a) + a V(x_i) ]
    """
    xm = np.roll(x, 1)
    kin = np.sum((x - xm) ** 2) / (4.0 * A)
    pot = A * np.sum(V(x))
    return float(kin + pot)


def grad_action(x: np.ndarray) -> np.ndarray:
    """
    Gradient of S w.r.t x_i, using periodic BC.

    From S_kin = sum_i (x_i - x_{i-1})^2/(4a):
      dS_kin/dx_i = (2 x_i - x_{i-1} - x_{i+1})/(2a)

    From S_pot = a sum_i V(x_i):
      dS_pot/dx_i = a V'(x_i)
    """
    xp = np.roll(x, -1)
    xm = np.roll(x, 1)
    gkin = (2.0 * x - xm - xp) / (2.0 * A)
    gpot = A * dVdx(x)
    return gkin + gpot


def gradient_cool(x0: np.ndarray, n_steps: int = 2000, eps0: float = 0.05) -> tuple[np.ndarray, list[float]]:
    """
    Deterministic cooling via gradient descent with a simple line search
    enforcing S_{new} <= S_{old}.
    """
    x = x0.copy()
    S_hist = [action(x)]
    eps = eps0

    for _ in range(n_steps):
        g = grad_action(x)

        S_old = S_hist[-1]
        # line search: shrink eps until action decreases
        step_ok = False
        eps_try = eps

        for _ls in range(20):
            x_try = x - eps_try * g
            S_new = action(x_try)
            if S_new <= S_old:
                x = x_try
                S_hist.append(S_new)
                step_ok = True
                # mildly increase step size if things are going well
                eps = min(eps_try * 1.05, 0.2)
                break
            eps_try *= 0.5

        if not step_ok:
            # cannot decrease action further at this resolution
            break

    return x, S_hist


def main():
    df = pd.read_csv(DATA / "quantum_path.csv")
    tau = df["tau"].to_numpy()
    x_raw = df["x"].to_numpy()

    # Choose steps automatically: enough to visibly smooth but not extreme.
    # For N=800 this is a good starting point.
    x_det, S_hist = gradient_cool(x_raw, n_steps=2500, eps0=0.05)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tau, x_raw, label="Monte Carlo", color="black", linewidth=1.0)
    ax.plot(tau, x_det, label="Deterministic cooling (Python)", color="mediumseagreen", linewidth=2.0)

    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$x(\tau)$")
    ax.set_title(r"Fig. 2 (det): Gradient-cooled path from raw MC")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()

    fig.savefig(FIGS / "Figure_2det.png", dpi=250)
    fig.savefig(FIGS / "Figure_2det.pdf")

    print(f"[ok] saved {FIGS/'Figure_2det.png'} and .pdf")
    print(f"[info] action: raw={action(x_raw):.6g}, det={action(x_det):.6g}, steps={len(S_hist)-1}")

    # Optional: save action history for debugging/tuning
    out_hist = pd.DataFrame({"step": np.arange(len(S_hist)), "S": S_hist})
    out_hist.to_csv(DATA / "det_cooling_action_history.csv", index=False)
    print(f"[ok] saved {DATA/'det_cooling_action_history.csv'}")


if __name__ == "__main__":
    main()