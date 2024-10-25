## Introduction
The code in this repositry is approximating the the time harmonic scattering of an incident plane wave $u^{i}(x) = e^{i k x \cdot d}$, where $k$ is the wavenumber and $d$ is the directional normal of the incident wave for $x \in \mathbb{R}^{2}$ on two disjoint 1D screens denoted $\Gamma_{j}$, for j = 1, 2. The solution $u$ satisfies
$$  \Delta u + k^{2} u = 0, ~~ x\in D:=\mathbb{R^{2}}\setminus \overline{\Gamma} ,
\\
u = 0, ~~ x \in \Gamma. $$
Additionally, the scattered field $u^{s} : = u - u^{i}$ satisfies the Sommerfeld radiation condition. This can be reformulated as the following Boundary integral problem,
$$ 	u(x) := u^{i}(x) - \frac{i}{4}\int_{\Gamma} H_{0}^{(1)}(k \vert x - y \vert) \phi (x), ~~ x \in D, $$
where $\phi$ is the solution to
$$ \frac{i}{4} \int_{\Gamma} H_{0}^{(1)} (k \vert x - y \vert ) \phi(y) d y = u^{i}(x), ~~ x\in \Gamma, $$
and $H_{0}^{(1)}$ is the Hankel function of the first kind.


Three different methods are given in this repository to solve this boundary integral equation using a boundary element method. A standard approach using a piecewise constant basis function is given in the folder *poly_approx_space*, as well as an iterative approach which also uses a piecewise constant basis functions. For details of the method please see [1].

An iterative Hybrid Numerical Asymptotic Boundary Element Method is given in the folder *HNABEM_iterative*. This uses an approximation spaces that is enriched with oscillatory basis functions, resulting in an increase in computational efficiency as $k$ grows. 

This code repository was used to produce the plots used in [1]. 

## Requirements
This code heavily uses the code HNABEMLAB by Andrew Gibbs [2, 3]. This has been slighty adapted for this purpose and the forked version ojp1995/HNABEMLAB (https://github.com/ojp1995/HNABEMLAB) should be used.
This code was ran using MATLAB 2022a, there might be some functions that may or may not work with newer or older versions. 

## References:
[1] O. Phillips, A Hybrid Numerical Asymptotic Boundary Element Method for Scattering by Multiple Screens, PhD Thesis, University of Reading, (2024)

[2] A. Gibbs, https://github.com/AndrewGibbs/HNABEMLAB

[3] A. Gibbs, D Hewett, D. Huybrechs, E. Parolin Fast hybrid numerical-asymptotic boundary element methods for high frequency screen and aperture problems based on least-squares collocation, SN Partial Differ. Equ. Appl, 1 (2020): 1-26.

