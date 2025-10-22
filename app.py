
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

st.set_page_config(page_title="Heat & Fluid Flow Playground – 2D Conduction", layout="centered")

st.title("Heat & Fluid Flow Playground – 2D Conduction")
st.write(
    """
This mini-app solves steady **2D heat conduction** (Laplace equation) in a rectangular domain
with Dirichlet left/right walls and adiabatic top/bottom. It’s a teaching sandbox: adjust parameters
and see how the temperature field and wall heat flux change.

**Model:** \\(\\nabla^2 T = 0\\) with left wall \\(T=T_{hot}\\), right wall \\(T=T_{cold}\\), top/bottom adiabatic (\\(\\partial T / \\partial y = 0\\)).
"""
)

with st.sidebar:
    st.header("Parameters")
    nx = st.slider("Grid points in x (nx)", 20, 200, 80, step=2)
    ny = st.slider("Grid points in y (ny)", 20, 200, 80, step=2)
    Lx = st.number_input("Width Lx", value=1.0, min_value=0.1, max_value=10.0, step=0.1, format="%.2f")
    Ly = st.number_input("Height Ly", value=1.0, min_value=0.1, max_value=10.0, step=0.1, format="%.2f")
    Thot = st.number_input("Left wall temperature Thot [K]", value=310.0, step=1.0, format="%.1f")
    Tcold = st.number_input("Right wall temperature Tcold [K]", value=290.0, step=1.0, format="%.1f")
    tol = st.number_input("Convergence tolerance", value=1e-6, format="%.1e")
    max_iter = st.number_input("Max iterations", value=5000, step=100)
    scheme = st.selectbox("Iteration scheme", ["Jacobi", "Gauss-Seidel"])
    sor_omega = st.slider("SOR omega (1=no SOR)", 1.0, 1.95, 1.2, step=0.05)
    st.caption("Tip: For Gauss-Seidel, omega in ~1.1–1.6 often converges faster.")

@st.cache_data(show_spinner=False)
def solve_conduction(nx, ny, Lx, Ly, Thot, Tcold, tol, max_iter, scheme, omega):
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)
    T = np.ones((ny, nx), dtype=float) * (Thot + Tcold) / 2.0

    # Dirichlet on left/right
    T[:, 0] = Thot
    T[:, -1] = Tcold

    dx2, dy2 = dx*dx, dy*dy
    beta = dx2 / dy2  # for compactness

    start = perf_counter()
    it = 0
    if scheme == "Jacobi":
        while it < max_iter:
            Told = T.copy()
            # interior update (Jacobi)
            T[1:-1,1:-1] = ( (Told[1:-1,2:] + Told[1:-1,:-2]) + beta*(Told[2:,1:-1] + Told[:-2,1:-1]) ) / (2*(1.0 + beta))
            # Neumann at top/bottom: dT/dy = 0 -> copy interior neighbor
            T[0,:]  = T[1,:]
            T[-1,:] = T[-2,:]
            # Dirichlet left/right
            T[:, 0]  = Thot
            T[:, -1] = Tcold
            it += 1
            err = np.max(np.abs(T - Told))
            if err < tol:
                break
    else:  # Gauss-Seidel with optional SOR
        while it < max_iter:
            err = 0.0
            # sweep interior
            for j in range(1, ny-1):
                for i in range(1, nx-1):
                    T_old = T[j, i]
                    T_new = ((T[j, i+1] + T[j, i-1]) + beta*(T[j+1, i] + T[j-1, i])) / (2*(1.0 + beta))
                    # SOR relaxation
                    T[j, i] = (1 - omega)*T_old + omega*T_new
                    diff = abs(T[j, i] - T_old)
                    if diff > err:
                        err = diff
            # Neumann top/bottom
            T[0,:]  = T[1,:]
            T[-1,:] = T[-2,:]
            # Dirichlet left/right
            T[:, 0]  = Thot
            T[:, -1] = Tcold
            it += 1
            if err < tol:
                break

    elapsed = perf_counter() - start

    # Heat flux (k=1) central differences
    qx = np.zeros_like(T)
    qy = np.zeros_like(T)
    qx[:,1:-1] = -(T[:,2:] - T[:,:-2])/(2*dx)
    qy[1:-1,:] = -(T[2:,:] - T[:-2,:])/(2*dy)

    q_left  = float(np.mean(qx[:, 0]))
    q_right = float(np.mean(qx[:, -1]))

    return T, qx, qy, q_left, q_right, it, elapsed

with st.spinner("Solving..."):
    T, qx, qy, q_left, q_right, it, elapsed = solve_conduction(nx, ny, Lx, Ly, Thot, Tcold, tol, int(max_iter), scheme, sor_omega)

st.subheader("Temperature field")
fig = plt.figure()
plt.contourf(T, levels=40)
plt.xlabel("x-index")
plt.ylabel("y-index")
st.pyplot(fig, clear_figure=True)

st.subheader("Heat flux (arrows)")
skip = max(1, int(max(nx, ny) / 25))
Y, X = np.mgrid[0:T.shape[0], 0:T.shape[1]]
fig2 = plt.figure()
plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], qx[::skip, ::skip], -qy[::skip, ::skip])
plt.xlabel("x-index")
plt.ylabel("y-index")
st.pyplot(fig2, clear_figure=True)

st.subheader("Convergence & Metrics")
st.write(f"- Iterations: **{it}**  |  Time: **{elapsed:.3f} s**")
st.write(f"- Mean heat flux at left wall: **{q_left:.4f}**  |  right wall: **{q_right:.4f}**")

with st.expander("Math & Notes"):
    st.markdown(
        r"""
**Governing equation (steady conduction):**  
\[\nabla^2 T = \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0\]

**Discretization (uniform grid):**  
\[(1+\beta)T_{i,j} = \tfrac{1}{2}\left(T_{i+1,j}+T_{i-1,j} + \beta (T_{i,j+1}+T_{i,j-1})\right),\quad \beta=\frac{\Delta x^2}{\Delta y^2}\]

**Boundary conditions:**  
- Left: \(T = T_{hot}\), Right: \(T = T_{cold}\) (Dirichlet)  
- Top/Bottom: \(\partial T/\partial y = 0\) (adiabatic)

**Heat flux (k=1):** \(\mathbf{q} = -\nabla T\)
"""
    )

st.caption("© 2025 Heat & Fluid Flow Playground • Educational use • MIT License")
