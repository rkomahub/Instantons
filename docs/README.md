# Monte Carlo Simulations of Instantons in the Double Well Potential

This project implements Euclidean lattice Monte Carlo simulations of instantons in the quantum mechanical double well potential.

The implementation follows:

T. Schäfer, *Instantons and Monte Carlo Methods in Quantum Mechanics*.

---

## 1. Physical Model

We study the Hamiltonian

\[
H = \frac{p^2}{2m} + (x^2 - \eta^2)^2
\]

After rescaling \(2m = 1\), the discretized Euclidean action reads

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

```
.
├── CMakeLists.txt
├── LICENSE
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

## 4. Build Pipeline (CMake)

### 4.1 Configure

Generate build files inside `build/`:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

### 4.2 Compile

Compile all source files and link executables:

```
cmake --build build -j 8
```

This produces:

```
bin/montecarlo     # Main simulation executable
bin/test_*         # Test executables
```

### 4.3 Execute Simulation

```
./bin/montecarlo
```

All numerical outputs are written as CSV files in the `data/` directory.

---

## 5. Running the Tests (CTest)

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

## 6. Documentation (Doxygen)

Requires `doxygen` (and optionally `graphviz`) to be installed.

Generate documentation:

```
cmake --build build --target docs
```

Open:

```
docs/api/html/index.html
```

Using wslview:

```
wslview $(pwd)/docs/api/html/index.html
```

Documentation is automatically extracted from comments inside `src/*.hpp` and `src/*.cpp`.

---

## 7. Reproducing Figures

After running the C++ executable, figures can be generated using:

```
python plots/plot_figX.py
```

The scripts reproduce the main figures of the Schäfer lectures:

* Potential and spectrum
* Euclidean paths and cooling
* Correlators and log-derivatives
* Instanton density vs cooling
* Instanton liquid comparisons

---

## 8. Scientific Goals

This code allows one to:

* Compare Monte Carlo results with semiclassical predictions
* Extract the energy gap from correlators
* Measure instanton density via cooling
* Study instanton liquid approximations
* Compute non-Gaussian tunneling corrections

---

## 9. Rebuild Options

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

Author: Marco Benazzi
Project started: 2026