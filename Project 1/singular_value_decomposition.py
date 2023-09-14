import numpy
from scipy.linalg import eig

def SVD(A : numpy.ndarray):
    """
    Single Value Decomposition Function. Performs SVD on a matrix A and returns
    1. Each matrix of the SVD decomposition
    2. The matrix condition number using its singular/eigenvalues
    3. The matrix inverse using the SVD decomposition
    """
    # first get eigenvalues of A*A^T and corresponding eigenvectors
    eigenvalues, _, V = eig(A*numpy.matrix.transpose(A)) # NOTE: matirx A is invertible if eigenvalues are all non-zero

    # converting eigenvalues to a diagonal matrix
    eigenvalues = numpy.diag(eigenvalues)
    
    if numpy.any(eigenvalues == 0):
        # matrix is not invertible
        raise Exception("Matrix is not invertible")

    # then calculate the singular values
    singular_values = numpy.sqrt(eigenvalues)

    # now calculate the U matrix using the formula U = A*V*sigma^-1
    U = A*V*numpy.linalg.inv(singular_values)

    # next, calculate condition number using singular values/eigenvalues
    # condition_number = 

    # finally, calculate the inverse using the SVD decomposition


    # return results
    return {
        "U": U,
        "sigma": singular_values,
        "V": V,
        "condition_number": numpy.linalg.cond(A),
        "inverse": numpy.linalg.inv(A)
    }