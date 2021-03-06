# Minimising computational cost for global search

## Overview
The minimisation of cost functions is crucial in various optimisation fields. However, identifying their global minimum remains challenging owing to the huge computational cost incurred. This work analytically expresses the computational cost to identify an approximate global minimum for a class of cost functions defined under a high-dimensional discrete state space. Then, we derive an optimal global search scheme that minimises the computational cost.

This project involves MATLAB scripts.

<br>

This is a supplementary material of

"Quadratic speedup of global search using a biased crossover of two good solutions"

Takuya Isomura

Preprint at https://arxiv.org/abs/2111.07680

<br>

Copyright (C) 2021 Takuya Isomura

(RIKEN Center for Brain Science)

<br>

2021-10-26


## System Requirements
This package requires only a standard computer with enough RAM to support the in-memory operations.

Software: MATLAB

RAM: 16+ GB

<br>

The package has been tested on the following system:

iMac Pro (macOS Catalina)

CPU: 2.3GHz 18 core Intel Xeon W

RAM: 128 GB

MATLAB R2020b

The runtimes below are generated using this setup.


## Demo
### fig1

The cost function and global minimum.

Relationship between the state dimensionality N and the global minimum L_GM.

Run 'fig1.m'. The simulation should take approximately 4 hours.

### fig2

Computational efficiency of gradient descent algorithms.

a, Schematic of basins and local minima in the N-dimensional state space.

b, Relationship between the state dimensionality N and the number of local minima N_LM.

Run 'fig2a.m' and 'fig2b.m'. The simulation should take approximately 1 day.

### fig3

Outcomes of selection and crossover algorithm.

Fourth-degree cost functions in a 200-dimensional binary state space (N = 200, K = 4) are employed, and M = 20000 states are sampled for each condition. The crossover rate is varied between 0 ≤ 𝛾 ≤ 0.5.

Run 'fig3.m'. The simulation should take approximately 38 hours.

### fig4

Application to the travelling salesman problem.

The number of cities is N = 500. The simulation is conducted with 50 different travelling salesman problems.

Run 'fig4.m'. The simulation should take approximately 28 hours.


## License
This project is covered under the **GNU General Public License v3.0**.
