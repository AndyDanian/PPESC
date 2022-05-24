
import numpy as np

def vector_to_matrix(n: int = None, vector: list = None, sym: str = None):
    """
    Transform one vector to symmetric or antisymmetric matrix

    Args:
        n(int) : dimension
        vector(list): vector array
        sym(str): symmetry of the matrix

    Return:
        matrix (list): array of 2D
    """

    if sym == "antisym":
        coef = -1.0
    else:
        coef = 1.0

    matrix = [[0 for i in range(n)] for j in range(n)]
    col = 0
    row = 0
    for x in vector:
        if sym == "square":
            matrix[row][col] = x
            col += 1
            if col == n:
                col = 0
                row += 1
        else:
            matrix[row][col + row] = x
            matrix[col + row][row] = coef * x
            col += 1
            if col == n - row:
                col = 0
                row += 1

    return matrix
