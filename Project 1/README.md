# Equilibirum Spring-Mass System Solver using Single Value Decomposition

## Table of Contents
1. [Overview](#overview)
2. [Single Value Decomposition](#single-value-decomposition)
3. [Usage](#usage)
    - [Fixed-Fixed System](#fixed-fixed-system)
    - [Fixed-Free System](#fixed-free-system)
    - [Free-Free System](#free-free-system)
4. [Discussion (Free-Free System)](#discussion-free-free-system)
5. [Dependencies](#dependencies)

## Overview
This is a Python program that solves a equilibrium spring-mass system using the [Single Value Decomposition (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition) method. It calculates and prints the displacements, elongations, and internal stresses in the system. For a more detailed explanation of the solver, read below.

## Equiliibrum Spring System Solver using SVD
More specifically, This program solves vertical spring mass systems like the one shown below (fixed-fixed system):
```
*wall*
__________
    | (spring)
    M (mass)
    | (spring)
    M (mass)
    | (spring)
__________    
*wall*
```

The program works by solving the linear system equations that models the spring-mass system. The system of equations is given by:

    f = Ku

Where f is the force vector (due to gravity), K is the stiffness matrix, and u is the displacement vector. The stiffness matrix is constructed using the spring constants and masses given by the user. Rearranging the terms allows us to solve for the displacement vector u:


    u = K^-1 * f
 
The solver then uses the SVD method to decompose the stiffness matrix into three matrices: `U`, `Σ`, and `V^T` and calculate `K^-1` using the equation below: 

    K^-1 = V * Σ^-1 * U^T

Then, the elongation vector is calculated by multiplying the transpose of V with the displacement vector. The internal stress vector is calculated by multiplying the stiffness matrix with the elongation vector. The solver then prints the displacement, elongation, and internal stress vectors.


## Usage
To run this program locally, you must have [Python3](https://www.python.org/downloads/) installed. Download [this github repository](https://github.com/Kelach/COE-352-Advanced-Scientific-Computation) and cd into this directory in your terminal.

Install the dependencies using the command below:

```
pip install -r requirements.txt
```

Then, run the program with the following command:
```
python equilibrium_solve.py
```
### Examples

#### Fixed-Fixed System (Two walls)
```
>> python equilibrium_solve.py
>> 
>> Enter your spring constants (N/m): 1 1
>> Enter your masses (kg): 1
>> 
>> SYSTEM SOLVED:
>> --------------------------------------------------
>> 
>> The displacements are: [4.905]
>> The elongations are: [ 4.905 -4.905]
>> The internal stresses are: [ 4.905 -4.905]
>> ==================================================
>> 
>> Condition number (L2): 1.0
>> Singular values: [2.]
>> Eigenvalues: [4.]
```
#### Fixed-Free System (One wall)
```
>> Enter your spring constants (N/m): 1 1 1 
>> Enter your masses (kg): 1 1 1
>> 
>> SYSTEM SOLVED:
>> --------------------------------------------------
>> 
>> ==================================================
>> 
>> Condition number: 16.39373162228439
>> Singular values: [3.2469796  1.55495813 0.19806226]
>> Eigenvalues: [10.54287655  2.41789479  0.03922866]
```

#### Free-Free System (No walls)
``` 
>> Enter your spring constants (N/m): 1
>> Enter your masses (kg): 1 1
>> ERROR - SYSTEM IS NOT SOLVABLE
```

## Discussion (Free-Free System)
You'll notice that the solver will return an error when the number of springs given is less than the number of masses. To explain why this is imagine the following:

Let's use the Free-Free System and above to construct our stiffness matrix. This will yield:

```
    --   --
    | 1 0 |
K = | 0 0 |
    --   --
```
If we wish to solve the system `f = Ku`. However, we see that our `K` matrix has a rank of one (there is only one linearly independent column/row) but requires a rank of 2. Therefore our stiffness matrix, in this example, is not full rank, and `K^-1` cannot be computed nor the spring-mass system be solved. In addition, we can also see that the `K` matrix is not full rank because its null space is not the zero vector. For example:

```
      --   --
      | 1 0 |
{0} = | 0 0 | * {0, 1}^T
      --   --
```
The null space of this matrix is vector span `[0 1]`. This means that the matrix `K` will map the vector `[0 1]` to the zero vector.

This aligns with our physical interpretations of the system too. A free-free spring-mass system implies that the masses connected to the springs are traversing through space freely. This means that the system is not constratined and the masses can move in any direction. Therefore, the position of the masses would be unknown.         

## Dependencies
This code relies on the following external libraries:

`numpy`: For numerical operations and array handling.

`scipy`: For eigenvalue calculations.
