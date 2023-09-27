import pytest
import numpy as np
def test_singular_value_decomposition():
    test_matrix1 = np.array([[5, -1],[5, 7]])
    test_matrix2 = np.array([[3, 3, 2],[2, 3, -2]])
    test_matrix3 = np.array([[3, 7, 6],[1, 2, 9],[2, 4, 8]])
    test_matrix4 = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
    test_matrix5 = np.array([[1, 2],[3, 4]])
    test_matrix6 = np.array([[1, 2, 3],[4, 5, 6]])
    test_matrix7 = np.array([[1, 2],[3, 4],[5, 6]])
    