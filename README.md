# Heat & Fluid Flow Playground – 2D Steady Heat Conduction

This project demonstrates the numerical solution of **steady-state two-dimensional heat conduction** in a rectangular domain using finite-difference methods.  
It provides an interactive visualization built in Python and Streamlit for educational and research purposes.

---

## 1. Physical Model

### 1.1 Governing Equation

For a homogeneous, isotropic medium with no internal heat generation, the steady conduction equation is:

\[
\nabla^2 T = 
\frac{\partial^2 T}{\partial x^2} +
\frac{\partial^2 T}{\partial y^2} = 0
\]

### 1.2 Domain and Boundary Conditions

- Geometry: rectangular domain of size \(L_x \times L_y\)
- Boundaries:
  - Left wall: \(T = T_\text{hot}\)
  - Right wall: \(T = T_\text{cold}\)
  - Top & bottom walls: adiabatic, \( \frac{\partial T}{\partial y} = 0 \)

The configuration represents **pure conduction** driven by a lateral temperature difference, with no convection or radiation.

---

## 2. Numerical Formulation

### 2.1 Finite-Difference Discretization

The Laplace equation is discretized on a uniform Cartesian grid \((i,j)\) as:

\[
T_{i,j}^{*} = 
\frac{
(T_{i+1,j} + T_{i-1,j}) + 
\beta \, (T_{i,j+1} + T_{i,j-1})
}{
2 (1 + \beta)
},
\quad
\beta = \frac{\Delta x^2}{\Delta y^2}
\]

This expression is derived from a second-order central-difference approximation of the spatial derivatives.

### 2.2 Iterative Solution

The system is solved iteratively using either:
- **Jacobi iteration**, or  
- **Gauss–Seidel iteration** (optionally with Successive Over-Relaxation, SOR)

The iteration continues until the maximum temperature difference between successive iterations satisfies:

\[
\max |T^{(n+1)} - T^{(n)}| < \varepsilon
\]

where \(\varepsilon\) is a user-defined convergence tolerance.

---

## 3. Post-Processing

### 3.1 Temperature Field

The steady-state temperature distribution \(T(x,y)\) is plotted as a filled contour map.

### 3.2 Heat-Flux Field

Heat flux is evaluated from Fourier’s law (assuming \(k = 1\)):

\[
\mathbf{q} = - \nabla T =
\left(
-\frac{\partial T}{\partial x},
-\frac{\partial T}{\partial y}
\right)
\]

and visualized as quiver arrows across the domain.

### 3.3 Wall Heat Fluxes

Average heat fluxes at the hot and cold walls are estimated from numerical gradients:

\[
\bar{q}_\text{left} =
\left\langle -\frac{\partial T}{\partial x} \right\rangle_{x=0}, 
\quad
\bar{q}_\text{right} =
\left\langle -\frac{\partial T}{\partial x} \right\rangle_{x=L_x}
\]

---

## 4. Educational Use and Interpretation

- Demonstrates the fundamentals of **steady conduction** without convection effects.  
- Helps visualize **isotherms** and **heat-flux vectors** under different mesh and relaxation parameters.  
- Forms a natural starting point for extending toward **transient conduction** or **natural convection**.

---

## 5. Validation

For the limiting case of one-dimensional conduction (small \(L_y\) or adiabatic top/bottom),  
the analytical solution reduces to a linear profile:

\[
T(x) = T_\text{hot} -
(T_\text{hot} - T_\text{cold}) \frac{x}{L_x}
\]

The mid-line numerical solution agrees closely with this analytical result.

---

## 6. Running the App

```bash
pip install -r requirements.txt
streamlit run app.py
