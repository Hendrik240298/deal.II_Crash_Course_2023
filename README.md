# deal.II Crash Course of the 10th deal.II Users and Developers Workshop

This repository contains the source code for the [deal.II](https://www.dealii.org) crash course of the [10th deal.II Users and Developers Workshop](https://www.dealii.org/workshop-2023/). 

The course is structured in two sessions. In the first session, the basic layout of a deal.II program will be covered by solving a stationary Poisson's problem. The second session will be dedicated to the solution of a time-dependent problem. More concrete, a linear Heat equation.

In both case, a code skeleton is provided which should be filled in a gap text manner. The sessions' folders can be found by their names. In addition, each session folder contains a solution folder which includes a working code solution.

## Session 1 
Session 1 aims to solve Poisson's problem and is oriented on [step-3](https://www.dealii.org/current/doxygen/deal.II/step_3.html) of the deal.II tutorial program. 

### Poisson's problem
Let $\Omega\subset\mathbb{R}^d$ and dimension $d=1,2,3$. Find $u:\bar{\Omega}\to\mathbb{R}$ such that

$$
\begin{aligned}
-\Delta u &= 1 &&\quad\text{in } \Omega\\
u &= 0 &&\quad\text{on } \partial\Omega
\end{aligned}
$$

Exemplary soltions for $d=1,2$ are sketched or shown below.
<div style="text-align:center">
    <img src="images/high_fidelity_1d.png" alt="Poisson 2D" height="175"/>
    <img src="images/poisson_2d.png" alt="Poisson 2D" height="175"/>
</div>

#### Exercise
The missing gaps are marked by `/* MISSING CODE /*` comments. The following steps must be done to complete the code:

1. Complete `Poisson::run()`.
   - This method is used as a central place to call the other methods of the class to solve the problem.
2. Complete `Poisson::make_grid()` by generating the mesh.
3. Complete `Poisson::assemble_system()`.
   1. Get the `FEValues<2> fe_values`.
   2. Reset the local cell's contributions.
   3. Assemble system matrix contributions.
   4. Assemble the rhs contribution. 
   5. Transfer the local elements to the global matrix.
   6. Transfer the local elements to the global rhs.
   7. Set Dirichlet boundary conditions.

#### Bonus
For further Interaction with the code, there are two bonus exercises that could be done to extend the code. These are marked in the code by `/* BONUS */`, 

1. Try out and replace the homogeneous Dirichlet BC by non-zero BC.
2. Investigate the solution on successively refined meshes.
