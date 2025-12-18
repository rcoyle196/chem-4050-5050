# CHEM 4050/5050: Computational Problem Solving (Fall 2025)

This repository contains homework submissions and project work focused on applying numerical methods, statistical analysis, and molecular simulations to chemical systems. All scripts were executed in a **Miniconda environment (v25.9.1)** utilizing **NumPy, SciPy, Matplotlib, and Pandas**.

---

## Homework 1: Introduction to Python
*Foundational scripts for data manipulation and physical chemistry constants.*

* **1-1: Ideal Gas Behavior** | Analyzing the relationship between volume and pressure in a closed system.
* **1-2: Database Integration** | Retrieving molecular geometry from the Computational Chemistry Comparison and Benchmark Database (CCCBDB).
* **1-Grad: Quantum Mechanics** | Solving the time-independent Schrödinger equation for a particle in a 1D box.

## Homework 2: Using Numerical Methods in Python
*Optimization and integration techniques applied to molecular potentials.*

* **2-1: Cluster Optimization** | Geometry optimization of Argon dimers and trimers using the Lennard-Jones potential.
* **2-2: Statistical Mechanics** | Numerical integration of virial coefficients.
* **2-Grad: Vibrational Spectroscopy** | Solving the Schrödinger equation for harmonic and anharmonic oscillators in one dimension.

## Homework 3: Statistics, Regressions, and Correlations
*Uncertainty analysis and linear modeling of thermodynamic properties.*

* **3-1: Trouton’s Rule** | Performing regression and uncertainty analysis on phase transition data.
* **3-Grad: Enthalpy Modeling** | Modeling the relationship between the enthalpy of vaporization and boiling points.

## Homework 4: Introduction to Molecular Simulation
*Thermodynamics and particle-level interactions.*

* **4-1: Computational Thermodynamics** | Computational analysis of thermodynamic work.
* **4-2: Lanthanide Properties** | Evaluating the thermodynamic properties of $Ce^{3+}$.
* **4-Grad: Dimer Dynamics** | Calculating the temperature of dimers using the Lennard-Jones potential.

## Homework 5: Monte Carlo Simulations
*Stochastic methods for quantum mechanical integration.*

* **5-1: Orbital Overlap** | Overlap integration of two Hydrogen $2p$ orbitals.
* **5-2: Kinetic Matrix Elements** | Diagonal and off-diagonal kinetic energy matrix elements between two Hydrogen $1s$ orbitals.

---

## Project 1: Competitive Adsorption
**Grand Canonical Monte Carlo (GCMC) Simulations**
Implementation of GCMC simulations to model the competitive adsorption of molecular species on surfaces, utilizing particle insertion, deletion, and displacement moves.

---

## Project 2: Polymer Chain Molecular Dynamics
**Folding & Unfolding Transitions**
This project simulates the conformational folding and unfolding of a 3D polymer chain.

### Parameter Configuration
**Standard Baseline:**
Unless specified otherwise, simulations used the following parameters:
* **Monomers ($N$):** 25
* **Spring Constant ($k$):** 5
* **Repulsive Epsilon ($\epsilon_{rep}$):** 0.25

**Parameter Sweeps Analyzed:**
* **Repulsion Strengths:** $\epsilon_{rep} = [0.25, 0.5, 0.75]$
* **Spring Constants:** $k = [5, 25, 50]$
* **Chain Lengths:** $N = [20, 25, 30, 50]$

### Directory Structure

* **End-to-End Distance:** Plots tracking head-to-tail separation, including derivative plots for transition temperature identification.
* **Potential Energy:** Total system potential energy across the thermal gradient, including derivative plots.
* **Radius of Gyration:** Plots quantifying the overall root-mean-square size and compactness of the chain.
* **Initial Velocities:** Histograms verifying the Maxwell-Boltzmann distribution at the simulation start.
* **Position Snapshots:** 3D visualizations of polymer conformations at $0.1T$, $0.5T$, and $1.0T$.

### Technical Notes
* **Units:** All measurements are reported in reduced units ($\sigma, \epsilon, \tau$).
* **Equilibration:** Data collection followed a 5,000-step equilibration period to ensure thermal stability.
