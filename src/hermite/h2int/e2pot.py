from lib2h import *


def e2pot(coord, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Electron repulsion integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        e2pot (ndarray): array 4d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    n: int = total_nprim

    #! Error with coord when only is one atom
    e2pot, count = i2e(
        np.asfortranarray(coord),
        np.array(lx),
        np.array(ly),
        np.array(lz),
        np.array(center),
        np.array(exp),
        dalton_normalization,
        n,
    )
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Electron Repulsion Atomic Integrals ({count})",
            delta_time=(time() - start),
        )

    return e2pot


if __name__ == "__main__":
    wf = wave_function(
        "../../tests/molden_file/LiH_sd.molden",
        scratch_path="/home1/scratch",
        job_folder="160922134451",
    )
    ## print wave function informatio
    print("\n Coordinates \n", wf.coordinates)
    print("\n # primtives \n ", wf.primitives_number)
    print("\n Exponents \n", wf.exponents)
    print("\n Centers \n", wf.primitives_centers)
    print("\n mlx \n ", wf.mlx)
    print("\n mly \n ", wf.mly)
    print("\n mlz \n ", wf.mlz)

    driver_time = drv_time()

    e2int = e2pot(
        coord=wf.coordinates,
        exp=wf.exponents,
        center=wf.primitives_centers,
        lx=wf.mlx,
        ly=wf.mly,
        lz=wf.mlz,
        output=11,
        dalton_normalization=False,
        driver_time=driver_time,
    )

    print("H2 STO2G : ", e2int)
    driver_time.printing()


#    print("\n\n",e2int)
# H2 ccpVTZ Integrals (160831) Time: 0:0:6.737
# H2 ccpVTZ Electron Repulsion Atomic Integrals (136288) Time: 0:0:6.335 con los 4 for
