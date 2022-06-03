from convert_array import *
import numpy as np

def print_triangle_matrix(integral: list = None, name: str = None, matriz_sym: str = None):
    """
    Print the triangule matrix

    Args:
        integral (array): array 2d with atomic integrals
        name (str): name of the integral
        matriz_sym (str): Matriz symmetric of atomic integrals
    """
    print("*" * 80)
    print("*** ", name.upper().center(70), " ***")
    print("*" * 80)

    ZERO = 1.0E-7

    size = len(integral[0][:])
    if size <= 5:
        chunks = size
        columns: int = chunks
    else:
        columns: int = 5
        chunks = int(size / 5)
        if size % 5 != 0:
            chunks += (size % 5) / (size % 5)

    count = 0
    while count < chunks:
        # Column number
        if count < chunks - 1:
            n = (count + 1) * columns
        else:
            n = size
        # Column index
        print(
            *[
                "    " + str(i + 1).center(14)
                if i == count * columns
                else str(i + 1).center(14)
                for i in range(count * columns, n)
            ],
            end="",
        )
        print()
        
        if matriz_sym == "square":
            initial_value: int = 0
        else:
            initial_value: int = count * columns

        for row in range(initial_value, size):
            # Row values and index
            if matriz_sym != "square":
                values = [integral[column][row]
                            for column in range(count * columns, n)
                            if row >= column
                            ]
            else:
                values = [integral[column][row]
                            for column in range(count * columns, n)
                            ]
            if np.linalg.norm(np.array(values)) > ZERO:
                print(
                    *[str(row + 1).center(4) 
                    + str("{:.6f}".format(value)).center(14)
                    if i == 0
                    else str("{:.6f}".format(value)).center(14)
                    for i, value in enumerate(values)],
                    end="",
                )
                print()

        if size < 5:
            count += chunks
        else:
            count += 1

def print_matriz_integrated(
    n: int = None, integrals: dict = None, symmetries: dict = None
):

    for integral_label in integrals.keys():
        print_triangle_matrix(
            vector_to_matrix(
                n,
                integrals[
                    integral_label
                ],
                symmetries[
                    integral_label
                ],
            ),
            integral_label,
            symmetries[
                integral_label
            ]
        )