# Continous Galerkin Finite Element Method (CG-FEM) for 1D Heat Equation in Python
## Table of Contents
1. [Introduction](#introduction)
2. [Mathematical Model](#mathematical-model)
3. [Weak Formulation](#weak-formulation)
4. [Update Equation Derivation](#update-equation-derivation)
    - [Forward Euler Method](#forward-euler-method)
    - [Backward Euler Method](#backward-euler-method)
6. [Results](#results)
    - [Forward Euler](#forward-euler)
    - [Backward Euler](#examples)
7. [Discussion](#discussion)
    - [Project Question 2](#project-question-2)
    - [Project Question 3](#project-question-3)

## Introduction
The finite element method (FEM) is a numerical technique for finding approximate solutions to boundary value problems for partial differential equations. It uses variational methods (the calculus of variations) to minimize an error function and produce a stable solution. Analogous to the idea that connecting many tiny straight lines can approximate a larger circle, FEM encompasses all the methods for connecting many simple element equations over many small subdomains, named finite elements, to approximate a more complex equation over a larger domain. FEM is a good choice for analyzing problems over complicated domains (like cars and oil pipelines), when the domain changes (as during a solid state reaction with a moving boundary), when the desired precision varies over the entire domain, or when the solution lacks smoothness. FEM allows for the solution of very large problems with reasonable computer resources; it is also the only approach for solving complex problems in, for example, 3D heat transfer or fluid dynamics as used in computer aided engineering, and nonlinear problems like implicit solvation as in computational chemistry. This code implements the continuous Galerkin finite element method for the 1D heat equation and applies Dirchelet Boundary Conditions.

- Scripts:
    - `main.py`: This script runs the CG-FEM approximation for the 1D heat equation using both the forward and backward Euler time-descretization methods. The results are compared to the analytical solution of the heat equation.
    `input.py`: This script contains the input parameters for the CG-FEM approximation.
    - `FEA.py`: This script contains the functions used to run the CG-FEM approximation for the 1D heat equation.
    - `experiment.py`: This script contains the functions used to run the experiments for the 1D heat equation.
    - `integrator.py`: This script contains the functions used to run the numerical integration for the 1D heat equation.

## Mathematical Model
The heat equation is a parabolic partial differential equation that describes the distribution of heat (or variation in temperature) in a given region over time. For a function u(x,t) that describes the temperature distribution over the rod, the heat equation this code solves is given by:
$$ u_t - u_{xx} = f(x,t) $$

where f(x, t) is given by:
$$ f(x,t) = (\pi^2 - 1)e^{-t}sin(\pi x) $$

This code is designed to solve this system with the followinkg initial and Dirchelet boundary conditions:
    
$$ u(x,0) = sin(\pi x) $$
$$ u(0,t) = 0 $$


The analytical solution to this system is given by:

$$ u(x,t) = e^{-t}sin(\pi x) $$


## Weak Formulation
The weak form of the heat equation is obtained by multiplying the equation by a test function $\psi$ and integrating over the domain $\Omega$. The test function is a function that is chosen to be sufficiently smooth so that the equation holds when it is multiplied by the test function. Multiplying both sides of the heat transfer equation gives:

$$
\frac{\partial u}{\partial t} \psi - \frac{\partial^2 u}{\partial x^2} \psi = f(x,t) \psi
$$

Integrating both sides of the equation over the domain `[0, 1]`:


$$
\int_{0}^{1} \frac{\partial u}{\partial t} \psi dx - \int_{0}^{1} \frac{\partial^2 u}{\partial x^2} \psi dx = \int_{0}^{1} f(x,t) \psi dx
$$

The second term on the left hand side of the equation can be integrated by parts:

$$
\int_{0}^{1} \frac{\partial u}{\partial t} \psi dx - \left.\frac{}{} \right|_{0}^{1} \frac{\partial u}{\partial x} \psi dx + \int_{0}^{1} \frac{\partial u}{\partial x} \frac{\partial \psi}{\partial x} dx = \int_{0}^{1} f(x,t) \psi dx
$$

After applying the boundary conditions the integrand of the second term on the left hand side of the equation goes to zero. Giving the final equation for the weak form:

$$
\int_{0}^{1} \frac{\partial u}{\partial t} \psi dx + \int_{0}^{1} \frac{\partial u}{\partial x} \frac{\partial \psi}{\partial x} dx = \int_{0}^{1} f(x,t) \psi dx
$$

## Update Equation Derivation
This code supports two the forward and backward euler time descritization methods. The update equation for each method is derived below.
### Forward Euler Method
Our weak form of the heat equation is given by:

$$
\int_{0}^{1} \frac{\partial u}{\partial t} \phi_i dx + \int_{0}^{1} \frac{\partial u}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f(x,t) \phi_i dx
$$
- Note:  $\psi_i = \phi_i$

Now, we approximate the time derivative of u(x, t) using the forward Euler method:

$$
\frac{\partial u}{\partial t} \approx \frac{u^{(n+1)} - u^{(n)}}{\Delta t}
$$

Substituting this approximation into our weak form of the heat equation gives:

$$
\int_{0}^{1} \frac{u^{(n+1)} - u^{(n)}}{\Delta t} \phi_i dx + \int_{0}^{1} \frac{\partial u^{n}}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f(x,t) \phi_i dx
$$

Now we apply Galerkin Expansion to the function `u` to get the following:

$$
\frac{1}{\Delta t} \int_{0}^{1} \sum_{j=1}^n u^{(n+1)}_j \phi_j \phi_i dx - \frac{1}{\Delta t} \int_{0}^{1} \sum_{j=1}^n u^{(n)}_j \phi_j \phi_i dx + \int_{0}^{1} \sum_{j=1}^n u^{(n)}_j \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f(x,t) \phi_i dx 
$$

Pulling out the summation gives:

$$
\frac{1}{\Delta t} \sum_{j=1}^n \int_{0}^{1} u^{(n+1)}_j \phi_j \phi_i dx - \frac{1}{\Delta t} \sum_{j=1}^n \int_{0}^{1} u^{(n)}_j \phi_j \phi_i dx + \sum_{j=1}^n \int_{0}^{1} u^{(n)}_j \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f(x,t) \phi_i dx
$$

We re-represent the integrals and summations with the following notation in matrix form:

- $M_{ij} = \int_{0}^{1} \phi_j \phi_i dx$

- $K_{ij} = \int_{0}^{1} \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx$

- $f^n_i = \int_{0}^{1} f(x,t) \phi_i dx$

- $u^{n} = \sum_{j=1}^n u^{(n)}_j$

- $u^{n+1} = \sum_{j=1}^n u^{(n+1)}_j$

Applyiing these subsitutions to the equation gives the following matrix form of the equation:

$$
\frac{1}{\Delta t} M u^{n+1} - \frac{1}{\Delta t} M u^{n} + K u^{n} = f^n
$$


Multiplying both sides by $\Delta t$:

$$
M u^{n+1} - M u^{n} + \Delta t K u^{n} = \Delta tf^n
$$

Moving the $M u^{n}$ and $\Delta t K u^{n}$ term to the right hand side:

$$
M u^{n+1} = M u^{n} - \Delta t K u^{n} + \Delta t f^n
$$

Multiplying both sides by $M^{-1}$:

$$
u^{n+1} = u^{n} - \Delta t M^{-1} K u^{n} + \Delta t M^{-1} f^n
$$

Finally, we can simplify the equation above to give the following update equation:

$$
u^{n+1} = (I - \Delta t M^{-1} K) u^{n} + \Delta t M^{-1} f^n
$$

### Backward Euler Method
To obtain the update equation for the backward Euler method, we follow the same steps as above, but approximate the time derivative of u(x, t) using the backward Euler method:

$$
\frac{\partial u}{\partial t} \approx \frac{u^{(n)} - u^{(n-1)}}{\Delta t}
$$

Substituting this approximation into our weak form of the heat equation gives:

$$
\frac{1}{\Delta t}\int_{0}^{1} u^{(n)} - u^{(n-1)} \phi_i dx + \int_{0}^{1} \frac{\partial u^{n}}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f^n \phi_i dx
$$

Shifting the time step by one to both sides of the equation gives:

$$
\frac{1}{\Delta t}\int_{0}^{1} u^{(n+1)} - u^{(n)} \phi_i dx + \int_{0}^{1} \frac{\partial u^{n+1}}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f^{n+1} \phi_i dx
$$

Now we apply Galerkin Expansion to the function `u` to get the following:

$$
\frac{1}{\Delta t} \int_{0}^{1} \sum_{j=1}^n u^{(n+1)}_j \phi_j \phi_i dx - \frac{1}{\Delta t} \int_{0}^{1} \sum_{j=1}^n u^{(n)}_j \phi_j \phi_i dx + \int_{0}^{1} \sum_{j=1}^n u^{(n+1)}_j \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f^{n+1} \phi_i dx
$$

Pulling out the summation gives:

$$
\frac{1}{\Delta t} \sum_{j=1}^n \int_{0}^{1} u^{(n+1)}_j \phi_j \phi_i dx - \frac{1}{\Delta t} \sum_{j=1}^n \int_{0}^{1} u^{(n)}_j \phi_j \phi_i dx + \sum_{j=1}^n \int_{0}^{1} u^{(n+1)}_j \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx = \int_{0}^{1} f^{n+1} \phi_i dx
$$

Once agian, we re-represent the integrals and summations with the following notation in matrix form:

- $M_{ij} = \int_{0}^{1} \phi_j \phi_i dx$

- $K_{ij} = \int_{0}^{1} \frac{\partial \phi_j}{\partial x} \frac{\partial \phi_i}{\partial x} dx$

- $f_i^n = \int_{0}^{1} f(x,t) \phi_i dx$

- $u^{n} = \sum_{j=1}^n u^{(n)}_j$

- $u^{n+1} = \sum_{j=1}^n u^{(n+1)}_j$

Applyiing these subsitutions to the equation gives the following matrix form of the equation:

$$
\frac{1}{\Delta t} M u^{n+1} - \frac{1}{\Delta t} M u^{n} + K u^{n+1} = f^{n+1}
$$


Moving the $\frac{1}{\Delta t} M u^{n}$ term to the right hand side:

$$
u^{n+1}(\frac{1}{\Delta t} M + K) = f^{n+1} + \frac{1}{\Delta t} M u^{n}
$$

Finally, Multiplying both sides by $(\frac{1}{\Delta t} M + K)^{-1}$ gives us the following update equation for the backward Euler method:

$$
u^{n+1} = (\frac{1}{\Delta t} M + K)^{-1} f^{n+1} + (\frac{1}{\Delta t} M + K)^{-1} \frac{1}{\Delta t} M u^{n}
$$

## Results
The results are compared to the analytical solution of the heat equation. The analytical solution is given by:

The solution approximations were calculated using the following configurations:
```
Number of Nodes: 11
Delta Time: 1/551
Spatial Bounds: [0, 1]
Time Bounds: [0, 1]
Boundary Conditions: [0, 'NA', 0, 'NA']
Quadrature Points: 2
Time Discretization Method: forwardEuler (or backwardEuler)

```

Below are the results of the code for the 1D heat equation using both Forward Euler and Backward Euler time-descretization methods. 
    - Note: The solution approximations are almost exact for both time-descretization methods due to the small time step used. Stability of the forward Euler method is further explored in the discussion section.

### Forward Euler
![Figure 1 Plot](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation/blob/main/Project%202/figures/figure1.png)

### Backward Euler
![Figure 2 Plot](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation/blob/main/Project%202/figures/figure2.png)


## Discussion
Solve first by using a forward Euler time derivative discretization
with a time-step of Δ𝑡 = 1
551 . Plot the results at the final time. Increase the
time-step until you find the instability. What dt does this occur at? How does
the solution change as N decreases?

### Project Question 2
The stability of the forward Euler method is explored by increasing the time step until the solution becomes unstable. The following plot demonstrates how our CG-FEM approximation becomes unstable as the time step increases for the forward Euler method:

![Figure 3](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation/blob/main/Project%202/figures/figure3.png)

From experimentation, I found that the smallest time step at which the approximation becomes unstable is 1/539. The following plot demonstrates how the solution changes as N decreases:

![Figure 5](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation/blob/main/Project%202/figures/figure5.png)

The plot shows that the solution becomes less precise as N decreases. This is because the solution is approximated using a smaller number of nodes, which means that the solution is approximated using a smaller number of basis functions leading to a lower fidelity approximation.

### Project Question 3
What happens as the time-step is equal to or greater than
the spatial step size? Explain why.
The following plot demonstrates how the solution changes as the time step is equal to or greater than the spatial step size for the backward Euler method:

![Figure 4](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation/blob/main/Project%202/figures/figure4.png)

Despite not diverging, the solution approximated using the backward Euler method becomes less accurate as the time step increases. This is because the backward Euler method is only first order accurate in time, meaning that the error in the solution will increase linearly with the time step size.