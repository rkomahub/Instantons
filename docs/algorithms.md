# Numerical Algorithms

This document describes the numerical procedures used in the instanton simulation project.

---

## 1. Lattice Discretization

Euclidean time is discretized as

\[
\tau_i = i a, \quad i = 0, \dots, N-1
\]

with periodic boundary conditions.

The discretized action is

\[
S = \sum_i
\left[
\frac{(x_i - x_{i-1})^2}{4a}
+ a (x_i^2 - \eta^2)^2
\right]
\]

---

## 2. Metropolis Algorithm

Each lattice site is updated sequentially.

Proposal:

\[
x_i \rightarrow x_i + \delta x,
\quad \delta x \sim \mathcal{N}(0, \sigma)
\]

Local action difference:

\[
\Delta S =
S_{\text{new}} - S_{\text{old}}
\]

Acceptance probability:

\[
P = \min(1, e^{-\Delta S})
\]

The proposal width is tuned to achieve ~50% acceptance.

---

## 3. Cooling Algorithm

Cooling modifies the update rule:

If

\[
\Delta S < 0
\]

then accept the update.

Else reject.

This removes ultraviolet fluctuations while preserving
large-scale semiclassical structures.

---

## 4. Correlator Computation

The two-point correlator is computed as

\[
C(\tau) =
\frac{1}{N}
\sum_i x_i x_{i+\tau}
\]

using periodic boundary conditions.

---

## 5. Ensemble Averaging

Given T independent configurations:

Mean:

\[
\bar{C}(\tau) =
\frac{1}{T}
\sum_{t=1}^T C_t(\tau)
\]

Unbiased sample variance:

\[
\mathrm{Var} =
\frac{1}{T-1}
\sum_t (C_t - \bar{C})^2
\]

Standard error:

\[
\mathrm{stderr} =
\sqrt{\mathrm{Var}/T}
\]

---

## 6. Instanton Counting

Instantons are identified by sign changes:

\[
x_i x_{i+1} < 0
\]

Total instanton density:

\[
n_{I+A} =
\frac{N_{I+A}}{\beta}
\]

---

## 7. Random Instanton Liquid Model (RILM)

Construct configuration:

\[
x(\tau) =
\eta
\left[
\sum_j Q_j
\tanh\left(\frac{\omega}{2}(\tau - \tau_j)\right)
- 1
\right]
\]

where positions τ_j are uniformly distributed.

---

## 8. Interacting Instanton Liquid Model (IILM)

Impose minimal separation:

\[
|\tau_i - \tau_j| \ge \tau_{\text{core}}
\]

This models short-range repulsion.

---

## 9. Non-Gaussian Correction (Adiabatic Switching)

Interpolate between Gaussian and full action:

\[
S_\alpha =
(1 - \alpha) S_{\text{gauss}}
+ \alpha S_{\text{full}}
\]

Compute integral numerically:

\[
\Delta S =
\int_0^1 d\alpha
\langle S_{\text{full}} - S_{\text{gauss}} \rangle_\alpha
\]

Integration method:
- Simpson rule (fine grid)
- Simpson rule (coarse grid)
- Richardson extrapolation:

\[
\Delta S \approx
\frac{16 I_{\text{fine}} - I_{\text{coarse}}}{15}
\]