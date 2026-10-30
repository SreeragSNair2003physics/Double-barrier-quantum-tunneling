
# Quantum Tunnelling Through a Double Barrier

**Numerical investigation of resonant quantum tunnelling using Crank–Nicolson time evolution and the transfer-matrix method**

## Overview

Quantum tunnelling allows a particle to penetrate a potential barrier even when its energy is below the barrier height. In a double-barrier system, repeated reflections between the barriers can produce resonant transmission at particular incident energies.

This project investigates one-dimensional electron tunnelling through a double rectangular potential barrier using two complementary numerical approaches:

1. **Time-dependent simulation:** Evolution of a Gaussian wave packet using the Crank–Nicolson finite-difference method.
2. **Stationary-state analysis:** Calculation of energy-dependent transmission and reflection using the transfer-matrix method.

The aim is to examine wave-packet dynamics, identify resonant transmission behaviour, and assess the numerical consistency and limitations of the calculations.

---

## 1. Physical Model

The system consists of two rectangular potential barriers separated by a central potential well.

### Simulation parameters

| Parameter | Value |
|---|---:|
| Particle | Electron |
| Barrier height | 0.8 eV |
| Barrier width | 0.5 nm |
| Central well width | 2.0 nm |
| Spatial domain | -30 nm to +30 nm |
| Spatial grid points | 5001 |
| Initial wave-packet centre | -4 nm |
| Initial wave-packet width parameter | 0.7 nm |
| Central wave number | $4\times10^9\ \mathrm{m}^{-1}$ |
| Time step | $5\times10^{-19}\ \mathrm{s}$ |
| Total simulated time | 30 fs |

The initial state is a normalized Gaussian wave packet with a complex phase corresponding to positive average momentum.

### Double-barrier potential

![Double-barrier potential and initial wave packet](Double%20barrier1.png)

---

## 2. Mathematical Framework


### 2.1 Time-dependent Schrödinger equation

The evolution of the wavefunction is governed by

$$
i\hbar\frac{\partial\psi(x,t)}{\partial t}=-\frac{\hbar^{2}}{2m}\frac{\partial^{2}\psi(x,t)}{\partial x^{2}}+V(x)\psi(x,t)
$$

Here, $\psi(x,t)$ is the complex wavefunction, $m$ is the electron mass, and $V(x)$ is the double-barrier potential.

The probability density is

$$
\rho(x,t)=\left|\psi(x,t)\right|^{2}
$$

For a normalized wavefunction, the total probability satisfies

$$
\int_{-\infty}^{\infty}\left|\psi(x,t)\right|^{2}\,dx=1
$$

### 2.2 Crank–Nicolson time evolution

The time-dependent simulation uses the Crank–Nicolson method, an implicit finite-difference scheme.

Writing the Schrödinger equation in terms of the Hamiltonian operator $H$, the time-stepping equation is

$$
\left(I+\frac{i\Delta t}{2\hbar}H\right)\psi^{n+1}=\left(I-\frac{i\Delta t}{2\hbar}H\right)\psi^{n}
$$

The spatial second derivative is discretized on a uniform grid:

$$
\frac{\partial^{2}\psi}{\partial x^{2}}\approx\frac{\psi_{j+1}-2\psi_j+\psi_{j-1}}{(\Delta x)^{2}}
$$

This produces a tridiagonal discrete Hamiltonian. The resulting linear system is solved using the **Thomas algorithm**, an efficient solver for tridiagonal matrices.

The simulation tracks the probability density as the incident wave packet interacts with the two barriers.

### 2.3 Transfer-matrix method

The stationary-state calculation solves the time-independent Schrödinger equation,

$$
-\frac{\hbar^{2}}{2m}\frac{d^{2}\psi(x)}{dx^{2}}+V(x)\psi(x)=E\psi(x)
$$

Within each constant-potential region, the wavefunction can be represented as a superposition of right-moving and left-moving components.

Continuity of the wavefunction and its spatial derivative at each interface is used to relate the amplitudes across the double-barrier structure.

The transfer-matrix calculation evaluates transmission and reflection probabilities over the energy interval **0.05–0.95 eV**.

For identical free-particle regions on both sides of a real, conservative potential, probability conservation requires

$$
T(E)+R(E)=1
$$

This provides an important validation condition for the stationary calculation.


## 3. Numerical Results

### 3.1 Wave-packet evolution

The Crank–Nicolson simulation shows the time-dependent interaction of the incident Gaussian wave packet with the double-barrier structure.

The calculation records the probability distribution in the left region, barrier regions, central well, and right region.

![Time-dependent wave-packet evolution](wavepacket_evolution.png)

The evolving probability density illustrates partial reflection, penetration into the central region, and transmission through the structure.

### 3.2 Resonant transmission

The stationary transfer-matrix calculation produces an energy-dependent transmission spectrum.

![Energy-dependent transmission spectrum](Double%20barrier2.png)

Pronounced transmission peaks are associated with resonant behaviour of the double-barrier system.

The locations and widths of these peaks depend on the barrier geometry and the numerical energy resolution.

### 3.3 Comparison of numerical approaches

The project also compares transmission results obtained from the stationary transfer-matrix method and the time-dependent Crank–Nicolson calculation.

![Comparison of transmission calculations](Final-comparison.png)

These approaches describe related but distinct quantities: stationary transmission is defined for an incident energy, while a finite wave packet contains a distribution of momenta and energies.

A quantitative comparison therefore requires careful treatment of the incident energy distribution and a well-defined procedure for extracting transmitted probability after scattering.

The comparison is included as an exploratory numerical assessment rather than a claim of exact agreement.

### 3.4 Boundary reflections

Finite computational domains can produce artificial reflections when propagating wave packets reach the domain boundaries.

The following figures investigate this effect.

**With boundary reflections:**

![Wave-packet evolution with boundary reflections](With%20reflections.png)

**Without boundary reflections:**

![Wave-packet evolution without boundary reflections](Without%20reflections.png)

These results illustrate why boundary treatment and simulation duration must be considered when interpreting time-dependent scattering calculations.

---

## 4. Numerical Validation

Several checks were performed to assess numerical reliability.

| Validation check | Purpose |
|---|---|
| Initial normalization | Verify that the starting wavefunction has unit total probability |
| Probability conservation | Monitor total probability during time evolution |
| Finite-value checks | Identify invalid numerical results in the transmission scan |
| $T+R$ conservation | Test the stationary scattering calculation |
| Spatial-grid sensitivity | Examine changes in observables with spatial resolution |
| Time-step sensitivity | Examine changes in observables with time-step size |

Numerical convergence tests were performed by varying spatial resolution and time step. The tested observables changed only slightly over the examined settings.

These checks support numerical stability for the quantities investigated, but they do not establish the absence of all discretization, boundary, or extraction errors.

---

## 5. Reproducing the Results

### Requirements

- Python
- NumPy
- Matplotlib

Clone the repository:

```bash
git clone https://github.com/SreeragSNair2003physics/Double-barrier-quantum-tunneling.git
```

Enter the repository directory:

```bash
cd Double-barrier-quantum-tunneling
```

Install the required packages:

```bash
python -m pip install numpy matplotlib
```

Run the main simulation:

```bash
python Double_barrier.py
```

The main script generates numerical diagnostics and plots associated with the double-barrier simulation.

The repository also contains additional calculations and figures used to investigate transmission, numerical validation, and boundary effects.

---

## 6. Limitations

The main limitations of the present study are:

1. **One-dimensional approximation:** The model considers a simplified one-dimensional scattering system rather than a complete three-dimensional physical device.
2. **Idealized barriers:** The potential consists of rectangular barriers with fixed parameters.
3. **Finite computational domain:** Boundary reflections can influence long-time wave-packet evolution.
4. **Finite numerical resolution:** The spatial grid, time step, and energy scan introduce numerical approximation errors.
5. **Transmission comparison:** Stationary-state and finite-wave-packet transmission require consistent definitions before quantitative agreement can be established.
6. **Resonance accuracy:** Precise resonance locations and widths require a sufficiently fine energy scan.

These limitations are important when distinguishing numerical observations from general physical conclusions.

---

## 7. Possible Extensions

Future work could include:

- Quantitative extraction of transmitted and reflected wave-packet probabilities after scattering.
- Energy-distribution-weighted comparisons between time-dependent and stationary transmission.
- More systematic grid and time-step convergence studies.
- Improved absorbing-boundary conditions.
- Resonance-width estimation and investigation of tunnelling lifetimes.

---

## 8. Tools and Techniques

**Programming:** Python, NumPy, Matplotlib

**Numerical methods:** Finite differences, Crank–Nicolson time evolution, Thomas algorithm, transfer matrices

**Physics:** Quantum tunnelling, wave-packet dynamics, resonant transmission, probability conservation

---

## Conclusion

This project develops complementary numerical approaches to quantum tunnelling through a double-barrier potential.

The time-dependent simulation provides a picture of wave-packet scattering, while the stationary transfer-matrix calculation reveals the energy dependence of transmission and the emergence of resonant peaks.

Together with numerical validation and boundary-reflection investigations, the project demonstrates the implementation and assessment of computational methods for a quantum-mechanical scattering problem.

**Project status:** Core numerical implementation and initial validation completed. Further quantitative comparison between the two transmission approaches remains a potential extension.
