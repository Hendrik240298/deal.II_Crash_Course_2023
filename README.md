# deal.II Crash Course of the 10th deal.II Users and Developers Workshop

This repository contains the source code for the [deal.II](https://www.dealii.org) crash course of the [10th deal.II Users and Developers Workshop](https://www.dealii.org/workshop-2023/).


## Session 1 
Session 1 aims to solve Poisson's problem and is oriented on [step-3](https://www.dealii.org/current/doxygen/deal.II/step_3.html) of the deal.II tutorial program. 

### Poisson's problem
Let $\Omega\subset\mathbb{R}^d$ and dimension $d=1,2,3$. Find $u:\bar{\Omega}\to\mathbb{R}$ such that
$$
\begin{align*}
    -\Delta u &= 1 \quad\text{in } \Omega\\
    u &= 0 \quad\text{on } \partial\Omega.
\end{align*}
$$
Exemplary soltions for $d=1,2$ are sketched or shown below.
<div style="text-align">
    <img src="images/high_fidelity_1d.png" alt="Poisson 2D" height="125"/>
    <img src="images/poisson_2d.png" alt="Poisson 2D" height="140"/>
</div>


