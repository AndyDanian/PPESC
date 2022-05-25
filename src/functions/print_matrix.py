from convert_array import *

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

    size = len(integral[0][:])
    chunks = int(size / 5)
    chunks += (size % 5) / (size % 5)
    count = 0
    while count < chunks:
        # Column number
        if count < chunks - 1:
            n = (count + 1) * 5
        else:
            n = size
        # Column index
        print(
            *[
                "    " + str(i + 1).center(14)
                if i == count * 5
                else str(i + 1).center(14)
                for i in range(count * 5, n)
            ],
            end="",
        )
        print()
        if matriz_sym == "square":
            for row in range(0, size):
                print(
                    *[
                        str(row + 1).center(4)
                        + str("{:.6f}".format(integral[column][row])).center(14)
                        if i == 0
                        else str("{:.6f}".format(integral[column][row])).center(14)
                        for i, column in enumerate(range(count * 5, n))
                        #if row >= column
                    ],
                    end="",
                )
                print()
            count += 1
        else:
            for row in range(count * 5, size):
                print(
                    *[
                        str(row + 1).center(4)
                        + str("{:.6f}".format(integral[column][row])).center(14)
                        if i == 0
                        else str("{:.6f}".format(integral[column][row])).center(14)
                        for i, column in enumerate(range(count * 5, n))
                        if row >= column
                    ],
                    end="",
                )
                print()
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