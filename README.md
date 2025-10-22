
# Heat & Fluid Flow Playground – 2D Conduction (Streamlit)

A tiny teaching app that solves steady **2D heat conduction** in a rectangular domain (Dirichlet left/right, adiabatic top/bottom) using Jacobi or Gauss–Seidel (+optional SOR). Built for *Heat & Fluid Flow Playground*.


## Problem Summary – 2D Steady Heat Conduction

This page presents an interactive mini-application that solves and visualizes a **steady-state 2D heat-conduction** problem in a rectangular domain (Lx × Ly).  
Its goal is to build intuition for conduction and to show how boundary conditions affect the temperature field and heat flux distribution.

### Geometric Model
- Domain: `Lx × Ly` (default 1 × 1)
- Grid: `nx × ny` points (uniform spacing)

### Governing Equation
For homogeneous conduction with no internal heat generation:
\[
\nabla^2 T = \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0
\]

### Boundary Conditions
- Left wall: **Dirichlet** \( T = T_{hot} \)  
- Right wall: **Dirichlet** \( T = T_{cold} \)  
- Top and bottom walls: **Adiabatic** \( \frac{\partial T}{\partial y} = 0 \)

> This setup represents a pure conduction case — heat flows horizontally from the hot to the cold side, while the top and bottom are insulated.

### Numerical Method
- Spatial discretization: second-order **finite differences** on a uniform grid.  
- Example interior update formula:
\[
T_{i,j}^{*}=\frac{(T_{i+1,j}+T_{i-1,j}) + \beta (T_{i,j+1}+T_{i,j-1})}{2(1+\beta)}, \quad \beta=\frac{\Delta x^2}{\Delta y^2}
\]
- Solver: **Jacobi** or **Gauss–Seidel** (with optional SOR relaxation).  
- Convergence: iteration stops when the maximum temperature change is below the user-defined tolerance.

### Outputs
- **Temperature field:** filled-contour plot.  
- **Heat-flux vectors:** \(\mathbf{q} = -\nabla T\) shown as arrows.  
- **Summary metrics:** iteration count, CPU time, average wall heat flux.

### Educational Notes
- This example illustrates **pure conduction only** — no fluid motion.  
  Natural convection would require adding momentum equations with the Boussinesq approximation.  
- The **adiabatic** top/bottom walls help visualize horizontal isotherms and no-flux boundaries.  
- Grid refinement and the SOR parameter affect convergence and smoothness — good for sensitivity studies.

### Simple Validation
In the 1-D conduction limit (small Ly or adiabatic top/bottom), the centerline profile \(T(x)\) becomes linear.  
The computed temperature along the midline can be compared with this analytical solution.

---

## Running Locally
```bash
pip install -r requirements.txt
streamlit run app.py
