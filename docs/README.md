# Monte Carlo Simulations of Instantons in the Double Well Potential

This project implements Euclidean lattice Monte Carlo simulations of instantons in the quantum mechanical double well potential.

The implementation follows the framework of:

T. Schäfer, *Instantons and Monte Carlo Methods in Quantum Mechanics*.

---

## 1. Physical Model

We study the Hamiltonian

\[
H = \frac{p^2}{2m} + (x^2 - \eta^2)^2
\]

After rescaling \(2m = 1\), the discretized Euclidean action is

\[
S = \sum_i \left[
\frac{(x_i - x_{i-1})^2}{4a}
+ a (x_i^2 - \eta^2)^2
\right]
\]

with periodic boundary conditions.

The classical instanton solution is

\[
x_I(\tau) = \eta \tanh\left(\frac{\omega}{2}(\tau - \tau_0)\right),
\quad \omega = 4\eta
\]

with classical action

\[
S_0 = \frac{4\eta^3}{3}.
\]

---

## 2. Implemented Methods

The project includes:

- Metropolis sampling of the Euclidean path integral
- Cooling for instanton extraction
- Random Instanton Liquid Model (RILM)
- Interacting Instanton Liquid Model (IILM)
- Adiabatic switching for non-Gaussian corrections
- Two-point correlators and gap extraction
- Ensemble averaging with statistical errors

---

## 3. Repository Structure
.
├── CMakeLists.txt # Build configuration
├── LICENSE
├── data/ # Generated CSV output files
├── plots/ # Python scripts to reproduce figures
│ ├── plot_fig1.py
│ ├── plot_fig2_6.py
│ ├── plot_fig3.py
│ ├── plot_fig4.py
│ └── plot_fig7.py
└── src/
├── main.cpp # Entry point
├── analysis_driver.* # High-level analysis pipeline
├── lattice.* # Path storage and boundary conditions
├── potential.* # Double well potential
├── metropolis.* # Monte Carlo evolution
├── cooling_evolution.* # Cooling analysis
├── observables.* # Correlation functions
├── instanton.* # Zero-crossing counting
├── ensemble.* # Ensemble averaging and error analysis
├── heating.* # Gaussian fluctuations
├── rilm.* # Random instanton liquid model
├── iilm.* # Interacting instanton liquid model
├── qmidens.* # Adiabatic switching density correction
└── parameters.hpp # Simulation parameters
---

## 4. Build Instructions

```
mkdir build
cd build
cmake ..
make
./instantons
```

All numerical outputs are written as CSV files in the `data/` directory.

---

## 5. Reproducing Figures

After running the C++ executable, figures can be generated using:

```
python plots/plot_figX.py
```

The plotting scripts reproduce the main figures of the Schäfer lectures:
- Potential and spectrum
- Euclidean paths and cooling
- Correlators and log-derivatives
- Instanton density vs cooling
- Instanton liquid comparisons

---

## 6. Scientific Goals

This code allows one to:

- Compare full Monte Carlo results with semiclassical predictions
- Extract the energy gap from correlators
- Measure instanton density via cooling
- Study instanton liquid approximations
- Compute non-Gaussian corrections to tunneling rates

---

## 7. Future Improvements

Possible extensions include:

- Autocorrelation time analysis
- Improved statistical error treatment
- Continuum limit study
- Parallel Monte Carlo implementation
- Extension to field-theory analogues

---

Author: Marco Benazzi  
Project started: 2026