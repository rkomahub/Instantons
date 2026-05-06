# Algorithms and Numerical Methods

This document describes the numerical algorithms used in the project.

The goal is to reproduce instanton physics in the quantum mechanical double well potential using Euclidean lattice Monte Carlo methods.

---

## 1. Discretized Euclidean Action

The continuum Euclidean action is discretized on a lattice with \(N\) points and spacing \(a\).

The lattice action is

$$
S_E[x] =
\sum_i
\left[
\frac{(x_i - x_{i-1})^2}{4a}
+
a (x_i^2 - \eta^2)^2
\right]
$$

Periodic boundary conditions are imposed.

The previous lattice site is

$$
i_- =
\begin{cases}
N-1, & i=0 \\
i-1, & i>0
\end{cases}
$$

The next lattice site is

$$
i_+ =
\begin{cases}
0, & i=N-1 \\
i+1, & i<N-1
\end{cases}
$$

---

## 2. Lattice Configuration

A configuration is represented by the vector

$$
x = (x_0, x_1, \dots, x_{N-1})
$$

The class `Lattice` stores this path.

Two initializations are used:

- hot start: random values in \([-\eta,\eta]\)
- cold start: constant configuration \(x_i = \eta\)

---

## 3. Metropolis Sampling

The Metropolis algorithm samples configurations with probability

$$
P[x] \propto e^{-S_E[x]}
$$

At each site, a local proposal is generated.

The proposal is

$$
x_i' = x_i + \delta x
$$

where

$$
\delta x \sim \mathcal{N}(0,\sigma)
$$

The action difference is computed locally.

The local action contribution is

$$
S_i =
\frac{(x_i-x_{i-1})^2 + (x_{i+1}-x_i)^2}{4a}
+
aV(x_i)
$$

The proposal is accepted with probability

$$
P_{\text{accept}} =
\min(1,e^{-\Delta S})
$$

where

$$
\Delta S = S_i[x_i'] - S_i[x_i]
$$

One sweep updates all lattice sites once.

---

## 4. Cooling

Cooling uses the same local proposal mechanism as Metropolis sampling, but the acceptance rule is deterministic.

A proposal is accepted if it does not increase the action.

$$
\Delta S \le 0
$$

Cooling removes short-distance fluctuations and drives configurations toward nearby classical structures.

The code uses cooling to measure:

- number of instantons
- instanton density
- action per instanton

---

## 5. Instanton Counting

Instantons and anti-instantons are detected through zero crossings.

A crossing is counted when neighboring points have opposite sign.

$$
x_i x_{i+1} < 0
$$

The total number of crossings gives

$$
N_{I+A}
$$

The corresponding density is

$$
n_{I+A} = \frac{N_{I+A}}{\beta}
$$

where

$$
\beta = Na
$$

---

## 6. Action Measurement

The function `compute_action` evaluates the full lattice action.

The kinetic term is computed from nearest-neighbor differences.

The potential term is computed from

$$
V(x_i) = (x_i^2-\eta^2)^2
$$

The total action is the sum over all sites.

---

## 7. Correlation Functions

The code computes Euclidean correlators from a single path.

For power \(p\), the correlator is

$$
C_p(\tau) =
\frac{1}{N}
\sum_i
x_i^p x_{i+\tau}^p
$$

with periodic indexing.

The cases used most often are:

- \(p=1\): \( \langle x(0)x(\tau) \rangle \)
- \(p=2\): \( \langle x^2(0)x^2(\tau) \rangle \)
- \(p=3\): \( \langle x^3(0)x^3(\tau) \rangle \)

For \(p=2\), the connected correlator is

$$
C_2^{\text{conn}}(\tau)
=
C_2(\tau)
-
\langle x^2\rangle^2
$$

---

## 8. Ensemble Averaging

Several independent measurements are combined to estimate means and uncertainties.

For measurements \(O_t\), the sample mean is

$$
\bar{O}
=
\frac{1}{T}
\sum_{t=1}^{T}
O_t
$$

The sample variance is

$$
s^2 =
\frac{1}{T-1}
\sum_{t=1}^{T}
(O_t-\bar{O})^2
$$

The standard error is

$$
\sigma_{\bar{O}} =
\sqrt{\frac{s^2}{T}}
$$

The code uses these formulas for correlator errors and density errors.

---

## 9. Random Instanton Liquid Model

The Random Instanton Liquid Model generates semiclassical paths from instanton centers.

The centers are sampled in Euclidean time.

The path is built from hyperbolic tangent profiles.

The general form is

$$
x(\tau)
=
\eta
\left[
\sum_j Q_j
\tanh\left(
\frac{\omega}{2}(\tau-\tau_j)
\right)
-1
\right]
$$

where

$$
\omega = 4\eta
$$

The charges are

$$
Q_j = \pm 1
$$

The generated paths are then measured with the same correlator routines used for Monte Carlo configurations.

---

## 10. Heating

Heating adds Gaussian fluctuations to a semiclassical path.

The update is

$$
x_i \rightarrow x_i + \xi_i
$$

where

$$
\xi_i \sim \mathcal{N}(0,\sigma)
$$

This produces noisy configurations around instanton backgrounds.

The heated paths are used to compare:

- pure RILM configurations
- RILM configurations with quantum-like fluctuations
- full Monte Carlo configurations

---

## 11. Interacting Instanton Liquid Model

The IILM represents configurations through collective coordinates.

A configuration stores:

- instanton positions \( \tau_j \)
- charges \( Q_j \)

A hard-core condition is imposed.

$$
d(\tau_i,\tau_j) \ge \tau_{\text{core}}
$$

The periodic distance is

$$
d(\tau_i,\tau_j)
=
\min(|\tau_i-\tau_j|,\beta-|\tau_i-\tau_j|)
$$

A proposal moves one instanton position.

$$
\tau_j' =
\tau_j + \delta \tau
$$

The proposed position is wrapped into \([0,\beta)\).

The proposal is rejected immediately if it violates the core constraint.

---

## 12. Instanton–Anti-Instanton Interaction

The instanton–anti-instanton interaction is studied by constructing trial IA configurations.

For each separation, the action is computed.

The interaction is measured relative to two isolated instantons.

$$
\frac{S_{\text{int}}}{S_0}
=
\frac{S}{S_0}
-
2
$$

The code also uses gradient-flow-like relaxation to generate streamline configurations.

The resulting paths and action densities are written to CSV files for plotting.

---

## 13. Adiabatic Switching

Adiabatic switching computes free-energy differences by interpolating between two actions.

The interpolating action is

$$
S_\alpha =
(1-\alpha)S_{\text{gauss}}
+
\alpha S_{\text{full}}
$$

where

$$
0 \le \alpha \le 1
$$

For each value of \(\alpha\), the code samples the corresponding ensemble and measures

$$
\Delta S(\alpha)
=
\langle S_{\text{full}} - S_{\text{gauss}} \rangle_\alpha
$$

The integral is

$$
I =
\int_0^1 d\alpha \, \Delta S(\alpha)
$$

---

## 14. Simpson Integration

The integral over \(\alpha\) is evaluated using Simpson's rule.

For an odd number of points \(n_\alpha\), the spacing is

$$
h = \frac{1}{n_\alpha-1}
$$

The Simpson approximation is

$$
I
\approx
\frac{h}{3}
\left[
f_0
+
f_{n_\alpha-1}
+
4\sum_{\text{odd }i} f_i
+
2\sum_{\text{even }i\ne 0,n_\alpha-1} f_i
\right]
$$

where

$$
f_i = \Delta S(\alpha_i)
$$

---

## 15. Richardson Extrapolation

The code computes the switching integral on two grids:

- fine grid
- coarse grid

The improved estimate is

$$
I_{\text{Rich}}
=
\frac{16I_{\text{fine}} - I_{\text{coarse}}}{15}
$$

This reduces the leading discretization error of Simpson integration.

---

## 16. Non-Gaussian Corrected Density

The Gaussian instanton density is computed from the semiclassical formula.

The non-Gaussian correction modifies it by

$$
n_{\text{corrected}}
=
n_{\text{gauss}}
e^{-I_{\text{Rich}}}
$$

This quantity is used in the density comparison against cooling estimates.

---

## 17. CSV Output Strategy

The C++ executable writes raw numerical data to the `data/` directory.

Python scripts in `plots/` read these files and generate figures.

The code separates:

- simulation
- measurement
- CSV export
- plotting

This makes the numerical output reproducible and easy to inspect.

---

## 18. Reproducibility

The executable accepts a seed through the command line.

The option is

```bash
--seed <int>
````

If no seed is given, the default seed is used.

This ensures that runs can be repeated exactly.

---

## 19. Command Workflow

The general workflow is

```bash
./bin/montecarlo <command>
python3 plots/plot_figX.py
```

For example, Fig.7 data are generated with

```bash
./bin/montecarlo fig7
```

and plotted with

```bash
python3 plots/plot_fig7.py
```

---

## 20. Code Organization

The source code is organized into layers.

The `core/` layer contains the basic physics and Monte Carlo engine.

The `models/` layer contains semiclassical instanton models.

The `analysis/` layer contains figure-specific and experiment-specific workflows.

The `utils/` layer contains generic tools such as I/O, statistics, and periodic indexing.

```
```