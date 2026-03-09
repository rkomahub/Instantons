# Monte Carlo Instantons in the Double Well Potential

Numerical simulations of instantons in the quantum mechanical double well potential using **Euclidean lattice Monte Carlo**.

The implementation follows:

**T. Schäfer**
*Instantons and Monte Carlo Methods in Quantum Mechanics*

---

# Overview

This project studies tunneling and instanton physics using several complementary numerical approaches:

* Euclidean **Monte Carlo path integral sampling**
* **Cooling** to extract instanton content
* **Instanton liquid models** (RILM / IILM)
* **Adiabatic switching** for non-Gaussian corrections
* **Correlation functions** and energy gap extraction

The code aims to reproduce the figures and numerical experiments presented in Schäfer’s lecture notes.

---

# Physical Model

We consider the double well Hamiltonian

[
H = \frac{p^2}{2m} + (x^2-\eta^2)^2
]

After rescaling ((2m = 1)), the discretized Euclidean action becomes

[
S =
\sum_i
\left[
\frac{(x_i-x_{i-1})^2}{4a}
+
a(x_i^2-\eta^2)^2
\right]
]

with periodic boundary conditions.

The classical instanton solution is

[
x_I(\tau) =
\eta
\tanh
\left(
\frac{\omega}{2}(\tau-\tau_0)
\right),
\qquad
\omega = 4\eta
]

and the classical instanton action

[
S_0 = \frac{4\eta^3}{3}
]

---

# Implemented Methods

### Monte Carlo path integral

* Metropolis updates of lattice configurations
* Ensemble averaging
* Correlator measurements
* Energy gap extraction

### Instanton extraction

* Cooling algorithms
* Instanton density measurements
* Action per instanton

### Semiclassical approaches

* Random Instanton Liquid Model (RILM)
* Heated RILM
* Interacting Instanton Liquid Model (IILM)

### Beyond Gaussian semiclassics

* Adiabatic switching
* Thermodynamic integration
* Non-Gaussian corrections to instanton density

---

# Repository Structure

```
.
├── CMakeLists.txt
├── LICENSE
├── README.md
├── bin/              # compiled executables
├── build/            # CMake build directory
├── data/             # generated CSV data
├── figures/          # generated plots
├── plots/            # python plotting scripts
├── docs/             # documentation (Doxygen)
├── tests/            # unit tests
└── src/
    ├── main.cpp
    ├── analysis_driver.*
    ├── lattice.*
    ├── potential.*
    ├── metropolis.*
    ├── cooling_evolution.*
    ├── observables.*
    ├── instanton.*
    ├── ensemble.*
    ├── heating.*
    ├── rilm.*
    ├── iilm.*
    ├── qmidens.*
    └── parameters.hpp
```

---

# Build

### Configure

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

### Compile

```bash
cmake --build build -j
```

Generated executables:

```
bin/montecarlo
bin/test_*
```

---

# Command Line Interface

The main executable runs different numerical experiments.

```
./bin/montecarlo <command> [options]
```

Optional arguments

```
--seed <int>
--etas a,b,c
```

Example

```
./bin/montecarlo eta-scan --etas 1.0,1.2,1.4,1.6
```

---

# Main Commands

| Command             | Description                               |
| ------------------- | ----------------------------------------- |
| `basic`             | typical Euclidean path + cooled path      |
| `ensemble`          | correlators and probability distributions |
| `cooling-evolution` | instanton density vs cooling              |
| `fig7`              | averaged cooling analysis                 |
| `fig8`              | instanton density vs η                    |
| `fig9`              | switching paths                           |
| `fig11`             | Gaussian effective potential              |
| `rilm`              | random instanton liquid model             |
| `heated-rilm`       | heated RILM                               |
| `iilm`              | interacting instanton liquid              |
| `fig14`             | IA interaction and streamline             |
| `fig15`             | streamline paths                          |
| `fig16`             | IA separation histogram                   |
| `fig17`             | instanton positions                       |
| `qmidens`           | non-Gaussian corrections                  |
| `eta-scan`          | scan in η                                 |

---

# Generating Figures

Typical workflow

```
./bin/montecarlo <command>
python3 plots/plot_figX.py
```

Example

```
./bin/montecarlo fig7
python3 plots/plot_fig7.py
```

---

# Figure Reproduction Map

| Schäfer figure | Method                           |
| -------------- | -------------------------------- |
| Fig.1          | Schrödinger spectral calculation |
| Fig.2          | Monte Carlo path                 |
| Fig.3          | probability distribution         |
| Fig.4          | correlators                      |
| Fig.5          | adiabatic switching              |
| Fig.6          | cooled correlators               |
| Fig.7          | cooling evolution                |
| Fig.8          | instanton density vs η           |
| Fig.9          | switching paths                  |
| Fig.10         | RILM correlators                 |
| Fig.11         | Gaussian effective potential     |
| Fig.12         | RILM vs heated RILM              |
| Fig.13         | heated RILM correlators          |
| Fig.14–16      | IA interaction and streamline    |
| Fig.17         | instanton liquid positions       |

---

# Testing

Run all unit tests

```
ctest --test-dir build --output-on-failure
```

Tests verify

* action consistency
* Metropolis acceptance rule
* cooling monotonicity
* boundary conditions
* potential symmetry

---

# Documentation

Generate API documentation with

```
cmake --build build --target docs
```

Then open

```
docs/api/html/index.html
```

---

# Computational Layers

The project contains four conceptual layers.

### Spectral quantum mechanics

Direct diagonalization of the Hamiltonian.

### Full Monte Carlo path integral

Metropolis sampling of Euclidean paths.

### Semiclassical instanton extraction

Cooling and instanton density measurements.

### Instanton liquid models

Analytical multi-instanton configurations (RILM, IILM, streamline).

---

# Rebuild

Incremental rebuild

```
cmake --build build -j
```

Clean rebuild

```
rm -rf build
cmake -S . -B build
cmake --build build -j
```

---

# Author

Marco Benazzi
Theoretical Physics MSc Project
Started: 2025