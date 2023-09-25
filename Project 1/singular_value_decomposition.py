import numpy as np
from numpy import matmul, diag
from scipy.linalg import eigh, svd
from collections import defaultdict


def SVD(A : np.ndarray):
    """
    Single Value Decomposition Function. Performs SVD on a matrix A and returns
    1. Each matrix of the SVD decomposition
    2. The matrix condition number using its singular/eigenvalues
    3. The matrix inverse using the SVD decomposition
    Dimensions: U: MxM, sigma: MxN, V: NxN
    """
    # @TODO: test current function with test input
    M, N = A.shape
    tol = 1*10**-5 # tolerance for singular values to be considered 0
    results = defaultdict(dict) # return variable

    # first get U, V matricies and eigenvalues
    eigenvalues, V = eigh(matmul(A.T, A)) # eigh sorts eigenvalues in ascending order (only to for symmetric matrices)
    eigenvalues, U = eigh(matmul(A, A.T)) 
    
    # change eigenvalues to descending order and change V, U accordingly
    eigenvalues = np.flip(eigenvalues)
    V = np.flip(V, axis=1)  
    U = np.flip(U, axis=1) 
    sigma = np.zeros((M, N))

    isInvertible = all(abs(eigenvalues) > tol) # check if matrix is invertible

    # now calculate singular matrix
    for i in range(min(M, N)):
        if abs(eigenvalues[i]) > tol :
            # if value is lower than tol we treat it as 0
            sigma[i][i] = np.sqrt(eigenvalues[i])
    
    # calculate inverse of sigma for A inverse
    sigma_inverse = np.zeros((M, N))
    for i in range(min(M, N)):
        if sigma[i][i] != 0:
            sigma_inverse[i][i] = 1 / sigma[i][i]

    # calculate inverse of A if A is invertible
    if isInvertible:
        results["inverse"] = matmul(matmul(U.T, sigma_inverse), V)

    results["U"] = U
    results["sigma"] = sigma
    results["V"] = V
    try:
        if isInvertible:
            results["conditionNumber"] = diag(sigma).max() / diag(sigma).min()
        else:
            results["conditionNumber"] = diag(sigma).max() / diag(sigma).min()
    except ZeroDivisionError:
        # if dividion by 0, that means condition number should be extremely large
        results["conditionNumber"] = float("inf") # l2 norm condition number 
    return results

tall = np.array([[-3, 1],[6, -2], [6, -2]])
wide = np.array([[3, 2, 2],[2, 3, -2]])
square = np.array([[-3, 1, 2],[6, -2, 3], [6, -2, 1]])
tallMatricies = np.array([])
for i in range(5):
    tallMatrix = np.random.randint(1, 11, size=(4, 3))
    result = SVD(tallMatrix)
    # check if decomposition is correct
    print("matrix: ", tallMatrix)
    print("re-composed:", matmul(result["U"], matmul(result["sigma"], result["V"].T)))
# U, s, V = svd(A)