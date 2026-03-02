# Monte Carlo Simulations of Instantons in the Double Well Potential

This project implements Euclidean lattice Monte Carlo simulations of instantons in the quantum mechanical double well potential.

The implementation follows:

**T. Schäfer**, *Instantons and Monte Carlo Methods in Quantum Mechanics*.

---

## 1. Physical Model

We study the Hamiltonian

[
H = \frac{p^2}{2m} + (x^2 - \eta^2)^2
]

After rescaling (2m = 1), the discretized Euclidean action reads

[
S = \sum_i \left[
\frac{(x_i - x_{i-1})^2}{4a}

* a (x_i^2 - \eta^2)^2
  \right]
  ]

with periodic boundary conditions.

The classical instanton solution is

[
x_I(\tau) = \eta \tanh\left(\frac{\omega}{2}(\tau - \tau_0)\right),
\quad \omega = 4\eta
]

with classical action

[
S_0 = \frac{4\eta^3}{3}.
]

---

## 2. Implemented Methods

The project includes:

* Metropolis sampling of the Euclidean path integral
* Cooling for instanton extraction
* Random Instanton Liquid Model (RILM)
* Interacting Instanton Liquid Model (IILM)
* Adiabatic switching for non-Gaussian corrections
* Two-point correlators and gap extraction
* Ensemble averaging with statistical errors

---

## 3. Repository Structure

```
.
├── CMakeLists.txt
├── LICENSE
├── bin/                   # Compiled executables
├── data/                  # Generated CSV output
├── plots/                 # Python plotting scripts
│   ├── plot_fig1.py
│   ├── plot_fig2_6.py
│   ├── plot_fig3.py
│   ├── plot_fig4.py
│   └── plot_fig7.py
├── docs/                  # Generated documentation (Doxygen)
├── tests/                 # Unit tests (CTest)
└── src/
    ├── main.cpp           # Command dispatcher
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

## 4. Build Pipeline (CMake)

### 4.1 Configure

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

### 4.2 Compile

```
cmake --build build -j 8
```

This produces:

```
bin/montecarlo
bin/test_*
```

---

## 5. Command-Line Interface

The executable now requires a command specifying which analysis to run. Since we introduced a **command-dispatch CLI**:
* The executable requires a **command**.
* Each command corresponds to a specific numerical experiment.
* Runs are reproducible via `--seed`.
* Optional `--etas` for eta scans.

### Basic Usage

```
./bin/montecarlo <command> [options]
```

### Figure to Command

| Schäfer Fig.   | Physical Content                                 | CLI Command                                     | Generated Data                            | Python Script               |
| -------------- | ------------------------------------------------ | ----------------------------------------------- | ----------------------------------------- | --------------------------- |
| **Fig. 1(a)**  | Double-well potential + energy levels            | *(no C++ needed)*                               | — (pure Python HO diagonalization)        | `plot_fig1.py`              |
| **Fig. 1(b)**  | Spectrum vs η                                    | *(no C++ needed)*                               | — (pure Python HO diagonalization)        | `plot_fig1.py`              |
| **Fig. 2**     | Typical Euclidean path + cooled path             | `basic`                                         | `quantum_path.csv`<br>`cooled_path.csv`   | `plot_fig2_6.py`            |
| **Fig. 3**     | Probability distribution vs exact |ψ₀|²          | `ensemble`                                      | `ensemble_quantum*`                       | `plot_fig3.py`              |
| **Fig. 4**     | Correlators (uncooled) + log-derivative          | `ensemble`                                      | `ensemble_quantum*`                       | `plot_fig4.py`              |
| **Fig. 5**     | Free energy via adiabatic switching              | `qmidens`                                       | switching CSV (ΔF, errors)                | *(future plot script)*      |
| **Fig. 6**     | Correlators (cooled) + log-derivative            | `ensemble`                                      | `ensemble_cooled*`                        | `plot_fig2_6.py`            |
| **Fig. 7(a)**  | Instanton density vs cooling sweeps              | `cooling-evolution`                             | `instanton_density_vs_ncool.csv`          | `plot_fig7.py`              |
| **Fig. 7(b)**  | Action per instanton vs cooling                  | `cooling-evolution` *(requires action logging)* | *(needs S per sweep export)*              | *(future plot)*             |
| **Fig. 8**     | Instanton density vs η                           | `eta-scan`                                      | `cooling_eta_scan.csv`                    | *(future plot)*             |
| **Fig. 9**     | One-instanton sector during switching            | `qmidens` *(future extension)*                  | *(not yet exported)*                      | —                           |
| **Fig. 10**    | RILM correlators                                 | `rilm`                                          | `rilm_path.csv`<br>`rilm_correlator.csv`  | *(future dedicated script)* |
| **Fig. 11**    | Gaussian effective potential                     | *(not implemented)*                             | —                                         | —                           |
| **Fig. 12**    | RILM vs Heated RILM path                         | `rilm` + `heated-rilm`                          | `rilm_path.csv`<br>`rilm_heated_path.csv` | *(future dedicated script)* |
| **Fig. 13**    | Heated RILM correlators                          | `heated-rilm`                                   | `rilm_heated_correlator.csv`              | *(future script)*           |
| **Fig. 14–16** | IA interaction, streamline, separation histogram | *(not yet implemented)*                         | —                                         | —                           |

---

### Optional Arguments

```
--seed <int>
```

Sets the random number generator seed (default: 12345).
Ensures reproducibility.

```
--etas a,b,c
```

Only for `eta-scan`.
Example:

```
./bin/montecarlo eta-scan --etas 1.0,1.2,1.4,1.6
```

---

## 6. Example Workflows

### Reproduce Fig. 2 (Typical Euclidean Path)

```
./bin/montecarlo basic --seed 12345
python plots/plot_fig2_6.py
```

### Reproduce Instanton Density vs Cooling

```
./bin/montecarlo cooling-evolution --seed 12345
python plots/plot_fig7.py
```

### Perform an η Scan

```
./bin/montecarlo eta-scan --etas 1.0,1.2,1.4,1.6
```

---

## 7. Running the Tests (CTest)

Unit tests verify:

* Action locality
* Cooling monotonicity
* Acceptance rule consistency
* Boundary-condition invariance
* Potential symmetry

Run all tests:

```
ctest --test-dir build --output-on-failure
```

---

## 8. Documentation (Doxygen)

Generate documentation:

```
cmake --build build --target docs
```

Open:

```
docs/api/html/index.html
```

---

## 9. Reproducing Figures

After generating data with the appropriate command, figures can be produced using:

```
python plots/plot_figX.py
```

The plotting scripts reproduce the main figures from the Schäfer lectures:

* Potential and spectrum
* Euclidean paths and cooling
* Correlators and log-derivatives
* Instanton density vs cooling
* Instanton liquid comparisons
* Non-Gaussian corrections

---

## 10. Rebuild Options

Incremental rebuild:

```
cmake --build build -j 8
```

Full clean rebuild:

```
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j 8
```

---

## 11. Project computational layers

## 1. Spectral (Schrödinger) Layer

**Method:** Direct diagonalization of the Hamiltonian
**Numerics:** HO basis → matrix → eigenvalues/eigenvectors
**Path integral not used**

### Figures:

* **Fig. 1(a)** — Potential + first energy levels
* **Fig. 1(b)** — Spectrum vs η

These are purely spectral. No Monte Carlo.

---

## 2. Full Monte Carlo (Euclidean Path Integral) Layer

**Method:** Metropolis sampling of discretized Euclidean action
**Core object:** Typical lattice path ( x(\tau_i) )
**Includes cooling as a post-processing filter**

### Figures:

* **Fig. 2** — Typical Euclidean path + cooled path
* **Fig. 3** — Probability distribution (P(x)) vs (|\psi_0|^2)
* **Fig. 4** — Correlators (uncooled)
* **Fig. 6** — Correlators (cooled)
* **Fig. 7(a)** — Instanton density vs cooling sweeps
* **Fig. 7(b)** — Action per instanton vs cooling
* **Fig. 8** — Instanton density vs η

These all originate from *full Monte Carlo sampling*.

Cooling does **not** create a new computational layer — it is a deterministic smoothing procedure applied to MC configurations.

---

## 3. Beyond Gaussian Semiclassics

**Method:** Adiabatic switching / thermodynamic integration
**Purpose:** Compute non-Gaussian corrections to instanton density
**Mathematically:** Interpolate between Gaussian and full action

### Figure:

* **Fig. 5** — Free energy (F(T))

Although it uses MC sampling internally, it belongs to a different conceptual layer because:

* It computes a *free energy difference*
* It probes *beyond one-loop semiclassics*
* It does not rely on raw path observables

This is a separate theoretical layer.

---

## 4. Semiclassical Instanton Liquid Models

**Method:** Construct multi-instanton configurations analytically
**No path integral sampling**

### Figures:

* **Fig. 9** — One-instanton sector during switching
* **Fig. 10** — RILM correlators
* **Fig. 11** — Gaussian effective potential
* **Fig. 12** — RILM vs Heated RILM paths
* **Fig. 13** — Heated RILM correlators
* **Fig. 14–16** — IA interaction, streamline, separation histogram

These are semiclassical constructions:

* RILM
* IILM
* Streamline
* IA interaction potential

No Metropolis sampling of full path integral.

# Conceptual Hierarchy

If we order them by theoretical depth:

1. Spectral QM (exact Schrödinger)
2. Full Monte Carlo path integral
3. Cooling extraction of semiclassical objects
4. Semiclassical liquid approximations
5. Beyond-Gaussian corrections (switching)

---

Author: Marco Benazzi
Project started: 2025