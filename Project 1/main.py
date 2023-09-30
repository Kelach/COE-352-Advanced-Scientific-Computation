import numpy as np
from numpy import array, zeros, diag, eye, matmul
from SVD import SVD

# @TODO: fix matrix inverse signing issue + add docstrings w GPT 
def main():
    """
    """
    # collect user input
    springs = input("Enter your spring constants (N/m): ")
    masses = input("Enter your masses (kg): ")

    # process input
    masses = array([float(mass) for mass in masses.split(" ")])
    springs = array([float(spring) for spring in springs.split(" ")])

    # defines our principal dimensions
    M, N = len(springs), len(masses) 

    # constructing A matrix (elongations)
    A = zeros((M, N))
    A += A + eye(M, N, k=0) + -1*eye(M, N, k=-1) # adds elongations

    # constructing C matrix (constituent matrix)
    C = diag(springs)

    # constructing force vectors (forces)
    forces = masses*9.81

    # constructing stiffness matrix (K)
    K = matmul(A.T, matmul(C, A)) # stiffness matrix

    # solving system
    try:
        svdResults = SVD(K) # calls custom SVD function
        if svdResults.get("inverse") is None:
            raise Exception
    except:
        # if K matrix is not invertible, the system cannot be solved
        print(f"""ERROR - SYSTEM IS NOT SOLVABLE""")
        return
    displacements = matmul(svdResults.get("inverse"), forces) 
    elongations = matmul(A, displacements)
    internalStresses = matmul(C, elongations)

    print(f"""
SYSTEM SOLVED:
{"-"*50}

The displacements are: {displacements}
The elongations are: {elongations}
The internal stresses are: {internalStresses}
{"="*50}

Condition number (L2): {svdResults.get("conditionNumber")}
Singular values: {diag(svdResults.get("sigma"))}
Eigenvalues: {diag(svdResults.get("sigma")**2)}
""")

main()