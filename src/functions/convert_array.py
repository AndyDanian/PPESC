from typing import Union

import numpy as np


def vector_to_matrix(
    n: int, vector: Union[np.ndarray, list], sym: str = ""
) -> np.ndarray:
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
        coef: float = -1.0
    else:
        coef = 1.0

    matrix: np.ndarray = np.zeros((n, n), dtype=float)
    col: int = 0
    row: int = 0
    for x in vector:
        if sym == "square":
            matrix[row, col] = x
            row += 1
            if row == n:
                row = 0
                col += 1
        else:
            matrix[row, col + row] = x
            matrix[col + row, row] = coef * x
            col += 1
            if col == n - row:
                col = 0
                row += 1

    return matrix
