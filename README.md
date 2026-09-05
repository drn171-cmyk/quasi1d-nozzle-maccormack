# 🚀 Quasi-1D Compressible Nozzle Solver

<p align="center">
  <b>MacCormack Predictor–Corrector Finite-Difference Method</b><br>
  Inviscid, quasi-one-dimensional compressible flow through a converging–diverging nozzle
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Language-C%2B%2B17-00599C?style=for-the-badge&logo=c%2B%2B&logoColor=white" alt="C++17">
  <img src="https://img.shields.io/badge/Method-MacCormack-orange?style=for-the-badge" alt="MacCormack">
  <img src="https://img.shields.io/badge/Flow-Compressible-red?style=for-the-badge" alt="Compressible flow">
  <img src="https://img.shields.io/badge/Model-Quasi--1D-6f42c1?style=for-the-badge" alt="Quasi 1D">
</p>

---

## ✨ Overview

This project implements a **quasi-one-dimensional compressible-flow solver** for a converging–diverging nozzle using the **explicit MacCormack finite-difference scheme**.

The solver advances the conservative form of the quasi-1D Euler equations in time and includes:

- ⚡ **Predictor–corrector MacCormack discretization**
- 📐 Geometric source term from the varying nozzle area
- 🌀 Local CFL-based adaptive time step
- 🧨 Pressure-based artificial viscosity for shock stabilization
- 🚪 Inlet and subsonic-exit pressure boundary conditions
- 📉 Density-based convergence monitoring
- 📊 Automatic export of the final solution to `nozzleResults.csv`

The current test case imposes a back pressure of `pBack = 0.6784`, producing a **normal shock inside the divergent section**.

---

## 🧠 Governing Equations

The quasi-1D Euler equations are written in conservative form as

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \frac{\partial \mathbf{F}}{\partial x}
= \mathbf{S}
$$

with the conservative variables

$$
\mathbf{U}=
\begin{bmatrix}
\rho A \\
\rho u A \\
\rho E A
\end{bmatrix}
$$

and the source vector

$$
\mathbf{S}=
\begin{bmatrix}
0 \\
 p\frac{dA}{dx} \\
0
\end{bmatrix}.
$$

The nozzle geometry used in the code is

$$
A(x)=1+2.2(x-1.5)^2,
\qquad 0\leq x\leq3.
$$

For a calorically perfect gas,

$$
\gamma=1.4,
$$

and the non-dimensional equation of state used by the solver is

$$
 p=\rho T.
$$

The local Mach number is evaluated from

$$
M=\frac{u}{a},
\qquad a=\sqrt{T}.
$$

---

## 🔢 Spatial Discretization
 
The domain is divided into `N = 1000` grid points:
 
```math
\Delta x = \frac{L}{N-1}
```
 
with
 
```math
L = 3.
```
 
### Predictor step — Forward difference
 
For each interior node,
 
```math
\mathbf{U}_i^* = \mathbf{U}_i^n - \frac{\Delta t}{\Delta x}\left(\mathbf{F}_{i+1}^n - \mathbf{F}_i^n\right) + \Delta t\,\mathbf{S}_i^n.
```
 
The forward spatial difference is therefore
 
```math
\left(\frac{\partial \mathbf{F}}{\partial x}\right)_i \approx \frac{\mathbf{F}_{i+1} - \mathbf{F}_i}{\Delta x}.
```
 
---
 
## 🔄 Corrector Step — Backward Difference
 
Fluxes are recalculated using the predicted state and the corrector uses a backward difference:
 
```math
\mathbf{U}_i^{n+1} = \frac{1}{2}\left[\mathbf{U}_i^n + \mathbf{U}_i^* - \frac{\Delta t}{\Delta x}\left(\mathbf{F}_i^* - \mathbf{F}_{i-1}^*\right) + \Delta t\,\mathbf{S}_i^*\right].
```
 
Therefore,
 
```math
\left(\frac{\partial \mathbf{F}}{\partial x}\right)_i \approx \frac{\mathbf{F}_i - \mathbf{F}_{i-1}}{\Delta x}
```
 
is used during the corrector stage.
 

This **forward/backward pairing** is the defining spatial discretization of the explicit MacCormack scheme implemented here.
---

## 📐 Discretization of the Area Source Term

During the predictor step, the source derivative is evaluated using a centered difference:

$$
\left(\frac{dA}{dx}\right)_i
\approx
\frac{A_{i+1}-A_{i-1}}{2\Delta x}.
$$

Hence,

$$
S_{2,i}=p_i\frac{A_{i+1}-A_{i-1}}{2\Delta x}.
$$

For the corrector step, the implemented source discretization is backward:

$$
S_{2,i}^*
=p_i^*\frac{A_i-A_{i-1}}{\Delta x}.
$$

The continuity and energy equations have zero geometric source terms in this formulation.

---

## 🛡️ Artificial Viscosity / Shock Capturing

A basic second-order MacCormack scheme can develop strong odd–even oscillations around a shock. To stabilize the present **shock-containing case**, the code adds pressure-based artificial viscosity.

The sensor is

$$
\varepsilon_i
=C_x
\frac{\left|p_{i+1}-2p_i+p_{i-1}\right|}
{p_{i+1}+2p_i+p_{i-1}},
$$

with

$$
C_x=0.3.
$$

The predicted conservative variables are then corrected with a second-difference damping term:

$$
\mathbf{U}_i^*
\leftarrow
\mathbf{U}_i^*
+\varepsilon_i
\left(
\mathbf{U}_{i+1}^n-2\mathbf{U}_i^n+\mathbf{U}_{i-1}^n
\right).
$$

The same idea is applied again after the corrector, but using the predicted-pressure field and predicted conservative variables.

This localized dissipation is especially important near the normal shock, where the physical solution contains a discontinuity.

---

## ⏱️ Time-Step Selection

The time step is selected from a local CFL condition and the **minimum value over the whole grid** is used:

$$
\Delta t_i
=CFL\frac{\Delta x}{u_i+a_i},
$$

$$
\boxed{\Delta t=\min_i(\Delta t_i)}.
$$

For the current shock-containing case,

$$
CFL=0.1.
$$

Using the minimum over the entire nozzle is important because the most restrictive point is not necessarily the outlet when a shock and subsonic exit are present.

---

## 🚪 Boundary Conditions

### Inlet

The reservoir conditions are non-dimensionalized such that

$$
\rho_0=T_0=1.
$$

The inlet momentum is extrapolated from the interior:

$$
U_{2,0}=2U_{2,1}-U_{2,2}.
$$

The total-energy variable is then reconstructed from the inlet velocity.

### Outlet

For the present subsonic-exit case, the outlet pressure is imposed:

$$
 p_e=p_{back}=0.6784.
$$

Density and velocity are extrapolated:

$$
\rho_e=2\rho_{N-2}-\rho_{N-3},
$$

$$
 u_e=2u_{N-2}-u_{N-3}.
$$

Temperature follows from

$$
T_e=\frac{p_e}{\rho_e}.
$$

The conservative variables are finally reconstructed from $(\rho_e,u_e,T_e)$.

---

## 📈 Numerical Result

The supplied implementation was compiled and executed with the parameters above.

### Convergence

The solver reached the convergence criterion

$$
\max_i\left|\frac{\rho_i^{n+1}-\rho_i^n}{\Delta t}\right|<10^{-8}
$$

after approximately **298,132 iterations**.

The final reported value was

$$
\max |d\rho/dt| \approx 9.99\times10^{-9}.
$$

### Output field

The final solution shows the expected compressible-flow behavior: acceleration through the converging section, supersonic flow in the divergent section, and a strong normal-shock transition toward the subsonic exit.

![MacCormack nozzle solution](nozzle_output.png)

> **Representative shock location:** approximately `x ≈ 2.08` in the supplied numerical run. The plotted shock marker is based on the strongest local Mach-number gradient in the generated output.

---

## 🧾 Output File

The solver writes the final solution to:

```text
nozzleResults.csv
```

with the columns:

| Column | Description |
|---|---|
| `x` | Axial coordinate |
| `A` | Cross-sectional area |
| `rho` | Density |
| `u` | Velocity |
| `T` | Static temperature |
| `p` | Static pressure |
| `Mach` | Local Mach number |

---

## ⚙️ Main Numerical Parameters

| Parameter | Value |
|---|---:|
| Specific heat ratio, `γ` | 1.4 |
| Nozzle length, `L` | 3 |
| Grid points, `N` | 1000 |
| CFL number | 0.1 |
| Artificial viscosity, `Cx` | 0.3 |
| Back pressure, `pBack` | 0.6784 |
| Convergence tolerance | `1e-8` |
| Maximum iterations | 400,000 |
| Scheme | MacCormack Predictor–Corrector |

---

## 🧩 Numerical Algorithm

```text
┌───────────────────────────────┐
│ Initialize x, A, ρ, u, T     │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Compute local Δt from CFL     │
│ Δt = min[ CFL·Δx / (u+a) ]    │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Predictor                     │
│ Forward difference            │
│ U → U*                        │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Artificial viscosity          │
│ Shock stabilization           │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Recalculate F* and S*         │
│ from predicted state          │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Corrector                     │
│ Backward difference           │
│ U* → Uⁿ⁺¹                      │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Artificial viscosity          │
│ Corrector stage               │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Apply boundary conditions     │
│ Update ρ, u, T, p, Mach       │
└───────────────┬───────────────┘
                │
                ▼
        Converged? ── No ──┐
                │          │
               Yes         └──► repeat
                │
                ▼
┌───────────────────────────────┐
│ Export nozzleResults.csv      │
└───────────────────────────────┘
```

---

## ▶️ How to Run

Compile with a C++17-compatible compiler:

```bash
g++ -O2 -std=c++17 nozzle.cpp -o nozzle
```

Run:

```bash
./nozzle
```

The solver then generates:

```text
nozzleResults.csv
```

---

## 📚 Method Summary

| Component | Discretization / Model |
|---|---|
| Governing equations | Quasi-1D Euler equations |
| Predictor | 1st-order forward difference |
| Corrector | 1st-order backward difference |
| Predictor source term | 2nd-order central difference for `dA/dx` |
| Corrector source term | Backward difference for `dA/dx` |
| Time integration | MacCormack predictor–corrector |
| Time step | Global minimum CFL condition |
| Shock stabilization | Pressure-based artificial viscosity |
| Convergence | Maximum `|dρ/dt|` criterion |

---

## 👨‍💻 Project Purpose

The project is intended as a compact numerical demonstration of **compressible CFD fundamentals**, particularly the connection between:

**conservative formulation → finite-difference discretization → predictor/corrector time marching → shock stabilization → convergence → post-processing**.

It is deliberately written without a CFD library so that the numerical method remains transparent and easy to modify.

---

## 📌 Notes

The variables `Cp` and `R` are defined in the source for the air model but are not directly required by the implemented non-dimensional update equations.

The current case is specifically configured to demonstrate a **shock-containing, subsonic-exit solution**. Lower back pressure values can be used to recover a fully supersonic, shockless configuration as described in the source code comments.

---

## ⭐ If you found this useful

Feel free to fork the project, experiment with the CFL number, grid resolution, back pressure, and artificial-viscosity coefficient, and compare the resulting shock position and flow-field profiles.
