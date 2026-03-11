# Diffusion-Limited Aggregation Simulation in Cylindrical Geometry (C implementation).
(Also called Diffusion Limited Deposition)

High-performance simulation code for particle deposition and cluster growth on a one–dimensional substrate using random walkers. The implementation includes several computational optimizations and analysis tools for studying interface growth, fractal properties, and tree statistics of the resulting aggregates.

The code is written in **C** and designed for large-scale simulations with up to millions of particles.

---

# Overview

This program simulates the stochastic deposition of particles on a linear substrate using random walkers that diffuse until they attach to the growing structure.

The simulation produces a growing cluster and computes several statistical properties of the interface and the internal structure of the aggregate, including:

- Interface roughness
- Skewness and kurtosis of the height distribution
- RMS thickness of the deposit
- Fractal dimension using box-counting
- Tree structure statistics (mass vs height scaling)

The code is optimized for **large systems** using spatial grids and collision prediction techniques.

---

# Relation to Ballistic Deposition Model (BDM)

Part of the dynamics implemented in this code corresponds to the **Ballistic Deposition Model (BDM) with η = 1**.

In this limit:

- The sticking probability becomes uniform.
- Growth is governed by local geometric constraints.
- Particles attach upon first contact with the cluster.

This implementation allows studying the structural properties of deposits that arise from this growth regime and comparing them with diffusion-limited aggregation dynamics.

---

# Main Features

## Random walker deposition

Particles are launched from above the current cluster height and perform random walks until they collide with the aggregate.

## Periodic boundary conditions

The substrate is periodic in the horizontal direction.

## Spatial acceleration structures

Several grids are used to improve performance:

- **Coarse collision grid** for fast collision prediction
- **Vicinity grid** to estimate safe step sizes
- **Occupancy grid (Y_grid)** for structural analysis

These reduce the computational cost of neighbor detection and particle collisions.

## Tree detection using BFS

The internal structure of the aggregate is analyzed by detecting **trees connected to the substrate** using a **Breadth-First Search (BFS)** algorithm.

For each tree the code measures:

- Tree mass `N_tree`
- Maximum tree height `H_tree`

The scaling relation between tree mass and height is recorded.

## Fractal dimension analysis

The fractal dimension of the aggregate is computed using a **box-counting method** applied at selected simulation times.

Also using the log-log of tree heights and mass.

## Interface statistics

For each time step the following quantities are computed:

- Mean interface height
- Roughness `W²`
- Skewness
- Kurtosis

These quantities allow studying universality classes of surface growth.

## RMS thickness

The vertical spread of deposited particles is measured via the RMS thickness of the deposit.

---

# Simulation Parameters

The main parameters of the simulation are defined in the code:

```
L  = 131072     // Substrate width
Lm = 4096        // Maximum simulation height
T  = 80          // Time steps
N  = 10          // Number of samples

L up to 1048576 with +70gb of RAM computers

```

Other relevant parameters:

```
diameter = 2
D_MAX = 160
COLLISION_THRESHOLD = 4
```

Snapshots for structural analysis are taken at:

```
t = {10, 20, 30, 40, 50, 60, 70, 80}
```

---

# Output Data

The simulation produces several output files.

## Interface statistics

```
roughness_opt_L=..._run=...dat
skewness_opt_L=..._run=...dat
kurtosis_opt_L=..._run=...dat
mean_height_opt_L=..._run=...dat
```

Active region statistics:

```
roughness_active_opt_L=..._run=...dat
skewness_active_opt_L=..._run=...dat
kurtosis_active_opt_L=..._run=...dat
mean_height_active_opt_L=..._run=...dat
```

---

## RMS thickness

```
rms_growth_avg_L=..._run=...dat
```

---

## Fractal dimension

For each snapshot time:

```
fractal_dim_t=..._L=..._run=...dat
```

Contains box-counting data of `N(s)` vs box size `s`.

---

## Tree statistics

Raw tree data:

```
all_trees_raw_t=..._run=...dat
```

Averaged tree scaling:

```
tree_scaling_avg_t=..._L=..._run=...dat
```

Columns:

```
N_tree   H_tree_avg   Count
```

---

## Final particle configuration

```
particulas_final.dat
```

Contains coordinates of all particles in the final configuration.

---

# Compilation

Compile using `gcc`:

```
gcc -Wall -O3 -march=native -mcmodel=large simulation.c -lm -o simulation
```

Recommended flags for performance:

```
-Wall -O3 -march=native -ffast-math

-mcmodel=large (obligatory flag for memory use)
```

---

# Running the Simulation

Run the program using:

```
./simulation run_id
```

Example:

```
./simulation 1
```

The `run_id` is used to differentiate output files.

---

# Performance Considerations

The code includes several optimizations:

- spatial hashing grids
- collision prediction using quadratic solutions
- dynamic memory allocation for neighbor lists
- BFS traversal for tree identification
- periodic boundary handling

These allow simulations with **millions of particles** while keeping computational cost manageable.


# Author

Derlis Alfonso

Simulation code developed for research on stochastic growth processes and fractal aggregates.
