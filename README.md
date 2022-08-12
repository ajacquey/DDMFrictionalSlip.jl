# DDMFrictionalSlip

[![Build Status](https://github.com/ajacquey/DDMFrictionalSlip.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ajacquey/DDMFrictionalSlip.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ajacquey/DDMFrictionalSlip.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ajacquey/DDMFrictionalSlip.jl)

DDMFrictionalSlip is a julia implementation of the Displacement Discontinuity Method (DDM) for two-dimensional domains (one-dimensional fracture).
Main features:
* Choice of Piecewise Constant (PWC), Piecewise Linear Collocation (PWLC), and Piecewise Quadratic Collocation (PWQ) shape functions
* Multithreaded assembly and solve
* Flexible problem formulation
* Non-equally sized elements

This package discretize the quasi-static changes in stress (normal or shear)  $\tau$ expressed as a integral of the displacement discontinuity $\delta$:

$$
    \tau\left(x\right) = \tau_{0}  + \frac{\mu^{\prime}}{\pi} \int_{\Omega} \frac{1}{s - x} \frac{\partial \delta}{\partial s} ds.
$$

$\tau_{0}$ is here the initial stress and $\mu^{\prime}$ the effective shear modulus. The previous expression is discretized into:

$$
    \tau_{i} = \tau_{0} + E_{ij} : \delta_{j},
$$

where $E_{ij}$ is the elastic collocation matrix (dense matrix).

This package can be used to solve for systems of coupled equations which can be expressed in the following way:

$$
    R_{\tau} = \Delta \tau - f_{\tau}\left(\Delta \epsilon, \Delta \delta\right) = 0
$$
$$
    R_{\epsilon} = \Delta \epsilon - f_{\epsilon}\left(\Delta \epsilon, \Delta \delta\right) = 0
$$

where $\Delta \tau$ and $\Delta \sigma$ are the changes in shear and normal stress respectively, $\Delta \delta$ and $\Delta \epsilon$ the changes in slip and opening repectively, and the two functions $f_{\tau}$ and $f_{\epsilon}$ can be defined to account for applied stress, frictional constraints, and/or fluid pressure coupling.

The user needs to specify the two functions $f_{\tau}$ and $f_{\epsilon}$ together with their derivatives with respect to the displacement discontinuity variables to properly calculate the jacobian matrix of the problem.
Please see the test suite `test/` for examples of formulations.

Author: Dr. Antoine B. Jacquey