# Molecular Dynamics Simulation of a Lennard–Jones Fluid

Author: Ledicia Díaz Lago
Language: Fortran (core engine) + Python (analysis)

---------------------------------------------------------------------

OVERVIEW

This project implements a Molecular Dynamics (MD) simulation of a three-dimensional monoatomic fluid interacting via the Lennard–Jones potential in the NVE ensemble.

The repository is modular and includes:

- Generation of an initial configuration
- Single production simulations
- Storage of time series of thermodynamic observables
- Storage of particle positions, velocities and accelerations
- Statistical analysis (means, correlations)
- Post-processing and visualization in Python
- An experimental framework for launching and averaging multiple simulations

IMPORTANT:
Before running a production simulation, the program
md_initial_config_program.f90
must be executed to generate the initial configuration from which the MD simulation starts.

---------------------------------------------------------------------

REPOSITORY STRUCTURE

```
scripts/
    -- md_simulation_program.f90
    -- md_initial_config_program.f90
    -- md_one_run_analysis.py
    -- base/
        -- define_precision.f90
        -- md_types.f90
        -- random_numbers.f90
        -- read_input_files.f90
    -- physics/
        -- geometry_pbc.f90
        -- lj_potential_energy.f90
        -- verlet.f90
        -- thermodynamic_coefs.f90
    -- stats/
        -- stats_math.f90
        -- md_means.f90
        -- md_correlations.f90
    -- run_many_md_simulations/
        -- md_simulation.f90
        -- md_initial_config.f90
        -- run_simulation.f90
        Important:
            - Work in progress
            - md_simulation.f90 and md_initial_config.f90 inside this directory are modules,
            not standalone programs. They are intended to be called by run_simulation.f90
            to allow launching multiple simulations and averaging results.

inputs/
    -- input_simulation_parameters.txt
outputs/
build/
```


---------------------------------------------------------------------

PHYSICAL MODEL

Particles interact via the Lennard–Jones potential:

$$U(r) = 4 ( 1/r^12 - 1/r^6 )$$

All quantities are expressed in fully reduced Lennard–Jones units:

$$\sigma = 1 \quad \epsilon = 1 \quad m = 1 \quad k_B = 1$$

The natural Lennard–Jones time scale is:

$$\tau = \sigma \sqrt{(m / \epsilon)}$$

Since $\sigma = m = \epsilon = 1$ in reduced units,
time in the simulation is already expressed in units of $\tau$. All variables are therefore dimensionless.

---------------------------------------------------------------------

VELOCITY–VERLET INTEGRATION

The equations of motion are integrated using the velocity–Verlet algorithm.

This scheme is obtained from a second-order Taylor expansion of the position:

$$r(t + dt) = r(t) + v(t) dt + \frac{1}{2} a(t) dt^2 + O(dt^3)$$

$$r(t - dt) = r(t) - v(t) dt + \frac{1}{2} a(t) dt^2 - O(dt^3)$$

Using Newton’s equation $(a = F/m, \mathrm{with}\;  m = 1)$, and combining forward and backward Taylor expansions, one obtains the symmetric velocity–Verlet scheme:

1) Position update: $r(t + dt) = r(t) + v(t) dt + \frac{1}{2} a(t) dt^2$

2) Force recomputation at the updated positions

3) Velocity update: $v(t + dt) = v(t) + \frac{1}{2} [ a(t) + a(t + dt) ] dt$

Properties:

- Second-order accurate in time
- Time-reversible
- Symplectic
- Good long-term energy conservation in the NVE ensemble

---------------------------------------------------------------------

THERMODYNAMIC OBSERVABLES

During the simulation, the following time series are stored:

- Potential energy $U$
- Kinetic energy $E_{kin} = \frac{1}{2} mv^2$
- Total energy $E_{tot}= U+E_{kin}$
- Temperature $T = \frac{2E_{kin}}{f}$ where $f$ is the number of degrees of freedom.
- Pressure $P = \rho * T + \frac{1}{3V} \sum_{i<j} ( r_{ij} F_{ij})$ where $\rho = N / V$.
- Other thermodynamic coeficientes such as $C_v,\; c_v,\; C_p,\; c_p,\; \gamma,\; k_S,\; \kappa_s,\; K_T,\; \kappa_T,\; \alpha_P,\; \alpha_S,\; \alpha_T$

---------------------------------------------------------------------

STORED DATA

The simulation produces:

- Time series of observables
- Binary trajectory file containing:
    - Wrapped positions
    - Unwrapped positions
    - Velocities
    - Accelerations
- Final statistical summary
- Autocorrelation files

Arrays for positions, velocities and accelerations are dynamically allocated,
allowing simulations with different numbers of particles without recompilation.

---------------------------------------------------------------------

STATISTICAL ANALYSIS

The code computes:

- Running averages
- Variances
- Block averages
- Autocorrelation functions

The Python analysis script additionally computes:

- Mean squared displacement (MSD)
- Velocity autocorrelation function (VACF)
- Diffusion coefficient (Einstein and Green–Kubo methods)
- Radial distribution function g(r)