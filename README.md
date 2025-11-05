# Stochastic Ribosome Transport Simulation

This Julia-based framework simulates ribosome traffic on a 1D lattice using a Gillespie stochastic simulation algorithm. The model incorporates pausing dynamics, particle clustering, and transient transport behavior under periodic boundary conditions (PBC). It supports comparison with multiple analytical approximations (Wang, Greulich, transient models) and includes visualization tools for current, density, and error metrics.

## 🚀 Features

- **Gillespie Simulation Engine**  
  Models ribosome dynamics with active/paused states, jammed particles, and periodic boundary conditions.

- **Transient Analysis**  
  Generates time-resolved current and density profiles over the lattice.

- **Analytical Model Benchmarking**  
  Compares simulation results to theoretical predictions using relative/absolute error metrics.

- **Interactive Plotting**  
  Tools to visualize cluster distributions, error isolines, and particle class decomposition via `PlotlyJS.jl` and `Interact.jl`.

- **Parallel Execution**  
  Supports distributed computation for parameter sweeps across $k_-, k_+$, $\rho$, and $L$.

- **Modular Architecture**  
  Codebase split into clearly defined modules for core dynamics, clustering analysis, plotting, and simulation orchestration.

---

## 📁 Directory Structure

```bash
.
├── main.jl                     # Top-level entry script for simulations
├── src/
│   ├── GillespieCore.jl        # Core dynamics and event rules
│   ├── Transient_Basics.jl     # Transient simulation components
│   ├── Clustering.jl           # Cluster detection and statistics
│   ├── Utilities.jl            # Diagnostics and consistency checks
│   ├── GillespieSimulation.jl  # Full simulation loops
├── plotting/
│   ├── pbc_plotting.jl         # Contour/error plots
│   ├── interactive_plotting_helpers.jl  # Interactive dashboards
├── data/                       # Simulation output and JLD2 data
