from io import TextIOWrapper

import numpy as np


def print_triangle_matrix(
    f: TextIOWrapper, integral: np.ndarray, matriz_sym: str = ""
) -> None:
    """
    Print the triangule matrix

    Args:
        integral (array): array 2d with atomic integrals
        matriz_sym (str): Matriz symmetric of atomic integrals
    """
    ZERO: float = 1.0e-8

    if integral.shape[0] != integral.shape[1]:
        raise ValueError("Matriz isn't square")

    size: int = integral.shape[0]

    if size <= 5:
        chunks: int = size
        columns: int = chunks
    else:
        columns = 5
        chunks = int(size / 5)
        if size % 5 != 0:
            chunks += int((size % 5) / (size % 5))

    count: int = 0
    while count < chunks:
        # Column number
        if count < chunks - 1:
            n = (count + 1) * columns
        else:
            n = size
        # Column index
        line: str = ""
        for i in range(count * columns, n):
            if i == count * columns:
                line += "    " + str(i + 1).center(16)
            else:
                line += str(i + 1).center(16)
        f.write(line + "\n\n")
        # print(
        #     *[
        #         "    " + str(i + 1).center(16)
        #         if i == count * columns
        #         else str(i + 1).center(16)
        #         for i in range(count * columns, n)
        #     ],
        #     end="",
        # )
        # print()

        if matriz_sym == "square":
            initial_value: int = 0
        else:
            initial_value = count * columns

        for row in range(initial_value, size):
            # Row values and index
            if matriz_sym != "square":
                values = [
                    integral[row, column]
                    for column in range(count * columns, n)
                    if row >= column
                ]
            else:
                values = [integral[row, column] for column in range(count * columns, n)]
            line = ""
            if np.linalg.norm(np.array(values)) > ZERO:
                for i, value in enumerate(values):
                    if i == 0:
                        line += str(row + 1).center(4) + str(
                            "{:.8f}".format(value)
                        ).center(16)
                    else:
                        line += str("{:.8f}".format(value)).center(16)
                f.write(line + "\n")
                # print(
                #     *[str(row + 1).center(4)
                #     + str("{:.8f}".format(value)).center(16)
                #     if i == 0
                #     else str("{:.8f}".format(value)).center(16)
                #     for i, value in enumerate(values)],
                #     end="",
                # )
                # print()

        if size < 5:
            count += chunks
        else:
            count += 1


if __name__ == "__main__":
    with open("a.dat", "a") as f:
        print_triangle_matrix(f, np.array([[1, 2], [3, 1]]))
