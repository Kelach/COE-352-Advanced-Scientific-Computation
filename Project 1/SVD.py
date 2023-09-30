from math import sqrt
import numpy as np
from numpy import matmul, diag
from scipy.linalg import eigh, svd
from collections import defaultdict


def SVD(A : np.ndarray, tol = 1*10**-5):
    """
    Single Value Decomposition Function. Performs SVD on a matrix A and returns
    1. Each matrix of the SVD decomposition
    2. The matrix condition number using its singular/eigenvalues
    3. The matrix inverse using the SVD decomposition
    Dimensions: U: MxM, sigma: MxN, V: NxN
    """

    results = defaultdict(dict)
    M, N = A.shape

    # first get U, V matricies and eigenvalues
    eigenvalues, U = eigh(matmul(A, A.T))
    _, V = eigh(matmul(A.T, A))

    # sort eigenvalues and corresponding U and V eigenvectors in descending order
    eigenvalues = np.flip(eigenvalues)
    U = np.flip(U, axis=1)
    V = np.flip(V, axis=1)

    # calculates sigma matrix
    sigma = np.zeros((M, N))
    for i in range(M):
        if eigenvalues[i] > tol:
            sigma[i][i] = sqrt(eigenvalues[i])

    # calculates inverse of sigma matrix for A^-1
    sigmaInverse = np.zeros((M, N))
    for i in range(M):
        if sigma[i][i] != 0:
            sigmaInverse[i][i] = 1 / sigma[i][i]
    
    # make sure that the signs of the columns of U and V match correctly
    same_sign = np.sign(((matmul(A, V))[0] * (matmul(U, sigma)[0])))
    V = V * same_sign.reshape(1, -1)
    
    if all(eigenvalues > tol): 
    # only calculate A^-1 if all eigenvalues are greater than "zero"
        results["inverse"] = matmul(matmul(V, sigmaInverse), U.T)
    

    # return results
    conditionNumber = diag(sigma).max() / diag(sigma).min() if diag(sigma).min() != 0 else float("inf")
    results["conditionNumber"] = conditionNumber
    results["U"] = U
    results["sigma"] = sigma
    results["V"] = V
    return results


# used for testing

# tall1 = np.array([[-3, 1],[6, -2], [6, -2]])
# tall2 = np.array([[-3, 1],[6, -2], [6, 4]])
# wide = np.array([[3, 2, 2],[2, 3, -2]])
# square = np.array([[-3, 1, 2],[6, -2, 3], [6, -2, 1]])
# tallMatricies = np.array([])

# for i in range(5):
#     tallMatrix = np.random.randint(1, 11, size=(3, 4))
#     result = SVD(tallMatrix)
#     print("matrix: ", tallMatrix)
#     print("re-composed:", matmul(result["U"], matmul(result["sigma"], result["V"].T)))