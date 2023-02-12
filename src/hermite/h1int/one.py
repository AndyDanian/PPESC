from lib1h import *


def one(
    exp,
    output,
    driver_time,
):
    """
    Matriz fill with 1
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    one: list = [1.0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Matrix fill with one",
            delta_time=(time() - start),
        )

    return one

if __name__ == "__main__":
    # STO-2G
    s = one(
        exp=[
            6.1638450, 1.0971610, 0.2459160, 0.0623710,
            0.2459160, 0.2459160, 0.2459160,
            0.0623710, 0.0623710, 0.0623710,
            1.3097564, 0.2331360
        ],
        output=3,
        driver_time=None,
    )

    print("A_B X p : ", s, "\n", len(s))
