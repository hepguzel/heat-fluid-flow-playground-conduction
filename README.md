
# Heat & Fluid Flow Playground – 2D Conduction (Streamlit)

A tiny teaching app that solves steady **2D heat conduction** in a rectangular domain (Dirichlet left/right, adiabatic top/bottom) using Jacobi or Gauss–Seidel (+optional SOR). Built for *Heat & Fluid Flow Playground*.

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud
1. Push this folder to a public GitHub repo (e.g., `hepguzel/heat-fluid-flow-playground-conduction`).
2. Go to https://share.streamlit.io , connect your GitHub, select the repo, **main file: `app.py`**.
3. Click **Deploy**. You’ll get a public URL.

## Roadmap
- Boundary options (all-Dirichlet, mixed BCs).
- 1D analytical comparison for validation (pure conduction).
- Non-uniform conductivity k(x,y).
- Separate page: correlation-based natural convection (Nu=f(Ra,Pr)).
- Separate page: OpenFOAM results visualizer (upload CSV/VTK → contour + quiver).

## License
MIT
