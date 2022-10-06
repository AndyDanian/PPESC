from libint import *
from constants_cto_gto import *

"""
Author: Mgs Andy Zapata
"""


def cto_gto_h1(Mxyz: np.ndarray, TP_A: np.ndarray) -> np.ndarray:
    """
    Convert cartesian to spherical integrals

    Arg:
        Mxyz [array, float]: Cartesian integrals
        TP_A [array, string]: Symbol associated with main quantum number

    Return:
        Mrtp [array, float]: Spherical integrals
    """

    Nsph: int = 0  # number of spherical primitives
    Np: int = 0  # number of cartesian primitives

    for a in TP_A:
        if a == "s":
            Nsph += 1
            Np += 1
        if a == "p":
            Nsph += 3
            Np += 3
        if a == "d":
            Nsph += 5
            Np += 6
        if a == "f":
            Nsph += 7
            Np += 10
        if a == "g":
            Nsph += 9
            Np += 15
        if a == "h":
            Nsph += 11
            Np += 21
        if a == "i":
            Nsph += 13
            Np += 28

    Mtemp: np.ndarray = np.zeros((Np, Nsph))
    for irow in range(Np):
        icol: int = 0
        jcol: int = 0
        for iket in TP_A:
            if iket == "s":
                Mtemp[irow, icol] = Mxyz[irow, jcol]
                icol += 1
                jcol += 1
            if iket == "p":
                Mtemp[irow, icol] = Mxyz[irow, jcol]
                Mtemp[irow, icol + 1] = Mxyz[irow, jcol + 1]
                Mtemp[irow, icol + 2] = Mxyz[irow, jcol + 2]
                icol += 3
                jcol += 3
            if iket == "d":
                Mtemp[irow, icol] = Mxyz[irow, jcol + 1]
                Mtemp[irow, icol + 1] = Mxyz[irow, jcol + 4]
                Mtemp[irow, icol + 2] = (
                    -0.5 * Mxyz[irow, jcol]
                    - 0.5 * Mxyz[irow, jcol + 3]
                    + Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 3] = Mxyz[irow, jcol + 2]
                Mtemp[irow, icol + 4] = (
                    R3_2 * Mxyz[irow, jcol] - R3_2 * Mxyz[irow, jcol + 3]
                )
                icol += 5
                jcol += 6
            if iket == "f":
                Mtemp[irow, icol] = (
                    -R18_4 * Mxyz[irow, jcol + 1] - R10_4 * Mxyz[irow, jcol + 6]
                )
                Mtemp[irow, icol + 1] = Mxyz[irow, jcol + 4]
                Mtemp[irow, icol + 2] = (
                    R30_20 * Mxyz[irow, jcol + 1]
                    + R6_4 * Mxyz[irow, jcol + 6]
                    + R30_5 * Mxyz[irow, jcol + 8]
                )
                # ml : 0
                Mtemp[irow, icol + 3] = (
                    R5_03 * Mxyz[irow, jcol + 2]
                    + R5_03 * Mxyz[irow, jcol + 7]
                    + Mxyz[irow, jcol + 9]
                )
                Mtemp[irow, icol + 4] = (
                    R6_4 * Mxyz[irow, jcol]
                    + R30_20 * Mxyz[irow, jcol + 3]
                    + R30_5 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 5] = (
                    R3_2 * Mxyz[irow, jcol + 2] - R3_2 * Mxyz[irow, jcol + 7]
                )
                Mtemp[irow, icol + 6] = (
                    R10_4 * Mxyz[irow, jcol] + R18_4 * Mxyz[irow, jcol + 3]
                )
                icol += 7
                jcol += 10
            if iket == "g":
                Mtemp[irow, icol] = (
                    R5_2 * Mxyz[irow, jcol + 1] - R5_2 * Mxyz[irow, jcol + 6]
                )
                Mtemp[irow, icol + 1] = (
                    -R2_D3_4 * Mxyz[irow, jcol + 4] - R10_4 * Mxyz[irow, jcol + 11]
                )

                Mtemp[irow, icol + 2] = (
                    R35_14 * Mxyz[irow, jcol + 1]
                    + R35_14 * Mxyz[irow, jcol + 6]
                    + R7_D3_7 * Mxyz[irow, jcol + 8]
                )

                Mtemp[irow, icol + 3] = (
                    R14_D3_28 * Mxyz[irow, jcol + 4]
                    + R70_D3_28 * Mxyz[irow, jcol + 11]
                    + R70_7 * Mxyz[irow, jcol + 13]
                )
                # ml : 0
                Mtemp[irow, icol + 4] = (
                    D3_8 * Mxyz[irow, jcol]
                    + R105_D3_140 * Mxyz[irow, jcol + 3]
                    + R105_D3_35 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 4] += (
                    D3_8 * Mxyz[irow, jcol + 10]
                    + R105_D3_35 * Mxyz[irow, jcol + 12]
                    + Mxyz[irow, jcol + 14]
                )

                Mtemp[irow, icol + 5] = (
                    R70_D3_28 * Mxyz[irow, jcol + 2]
                    + R14_D3_28 * Mxyz[irow, jcol + 7]
                    + R70_7 * Mxyz[irow, jcol + 9]
                )
                Mtemp[irow, icol + 6] = (
                    R5_4 * Mxyz[irow, jcol]
                    + R21_D3_14 * Mxyz[irow, jcol + 5]
                    - R5_4 * Mxyz[irow, jcol + 10]
                )
                Mtemp[irow, icol + 6] -= R21_D3_14 * Mxyz[irow, jcol + 12]
                Mtemp[irow, icol + 7] = (
                    R10_4 * Mxyz[irow, jcol + 2] + R2_D3_4 * Mxyz[irow, jcol + 7]
                )
                Mtemp[irow, icol + 8] = (
                    R35_8 * Mxyz[irow, jcol]
                    + R3_D3_4 * Mxyz[irow, jcol + 3]
                    + R35_8 * Mxyz[irow, jcol + 10]
                )
                icol += 9
                jcol += 15
            if iket == "h":
                # ml: -5
                Mtemp[irow, icol] = (
                    H910 * Mxyz[irow, jcol + 1]
                    + H93 * Mxyz[irow, jcol + 6]
                    + H90 * Mxyz[irow, jcol + 15]
                )
                # ml: -4
                Mtemp[irow, icol + 1] = (
                    H39 * Mxyz[irow, jcol + 4] - H39 * Mxyz[irow, jcol + 11]
                )
                # ml: -3
                Mtemp[irow, icol + 2] = (
                    H50 * Mxyz[irow, jcol + 1]
                    - H53 * Mxyz[irow, jcol + 6]
                    - H512 * Mxyz[irow, jcol + 8]
                )
                Mtemp[irow, icol + 2] = (
                    Mtemp[irow, icol + 2]
                    - H50 * Mxyz[irow, jcol + 15]
                    - H55 * Mxyz[irow, jcol + 17]
                )
                # ml: -2
                Mtemp[irow, icol + 3] = (
                    H44 * Mxyz[irow, jcol + 4]
                    + H44 * Mxyz[irow, jcol + 11]
                    + H114 * Mxyz[irow, jcol + 13]
                )
                # ml: -1
                Mtemp[irow, icol + 4] = (
                    H110 * Mxyz[irow, jcol + 1]
                    + H13 * Mxyz[irow, jcol + 6]
                    + H112 * Mxyz[irow, jcol + 8]
                )
                Mtemp[irow, icol + 4] += (
                    H10 * Mxyz[irow, jcol + 15]
                    + H15 * Mxyz[irow, jcol + 17]
                    + H114 * Mxyz[irow, jcol + 19]
                )
                # ml: 0
                Mtemp[irow, icol + 5] = (
                    H02 * Mxyz[irow, jcol + 2]
                    + H07 * Mxyz[irow, jcol + 7]
                    + H09 * Mxyz[irow, jcol + 9]
                )
                Mtemp[irow, icol + 5] += (
                    H02 * Mxyz[irow, jcol + 16]
                    + H09 * Mxyz[irow, jcol + 18]
                    + Mxyz[irow, jcol + 20]
                )
                # ml: 1
                Mtemp[irow, icol + 6] = (
                    H10 * Mxyz[irow, jcol]
                    + H13 * Mxyz[irow, jcol + 3]
                    + H15 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 6] += (
                    H110 * Mxyz[irow, jcol + 10]
                    + H112 * Mxyz[irow, jcol + 12]
                    + H114 * Mxyz[irow, jcol + 14]
                )
                # ml: 2
                Mtemp[irow, icol + 7] = (
                    H32 * Mxyz[irow, jcol + 2]
                    + H39 * Mxyz[irow, jcol + 9]
                    - H32 * Mxyz[irow, jcol + 16]
                    - H39 * Mxyz[irow, jcol + 18]
                )
                # ml: 3
                Mtemp[irow, icol + 8] = (
                    H50 * Mxyz[irow, jcol]
                    + H53 * Mxyz[irow, jcol + 3]
                    + H55 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 8] += (
                    H512 * Mxyz[irow, jcol + 12] - H50 * Mxyz[irow, jcol + 10]
                )
                # ml: 4
                Mtemp[irow, icol + 9] = (
                    H72 * Mxyz[irow, jcol + 2]
                    + H77 * Mxyz[irow, jcol + 7]
                    + H72 * Mxyz[irow, jcol + 16]
                )
                # ml: 5
                Mtemp[irow, icol + 10] = (
                    H90 * Mxyz[irow, jcol]
                    + H93 * Mxyz[irow, jcol + 3]
                    + H910 * Mxyz[irow, jcol + 10]
                )  # C55
                icol += 11
                jcol += 21
            if iket == "i":
                # ml: -6
                Mtemp[irow, icol] = (
                    I12_1 * Mxyz[irow, jcol + 1]
                    + I12_6 * Mxyz[irow, jcol + 6]
                    + I12_15 * Mxyz[irow, jcol + 15]
                )
                # ml: -5
                Mtemp[irow, icol + 1] = (
                    I10_4 * Mxyz[irow, jcol + 4]
                    + I10_11 * Mxyz[irow, jcol + 11]
                    + I10_22 * Mxyz[irow, jcol + 22]
                )
                # ml: -4
                Mtemp[irow, icol + 2] = (
                    I8_1 * Mxyz[irow, jcol + 1]
                    + I8_8 * Mxyz[irow, jcol + 8]
                    + I8_15 * Mxyz[irow, jcol + 15]
                )
                Mtemp[irow, icol + 2] += I8_17 * Mxyz[irow, jcol + 17]
                # ml: -3
                Mtemp[irow, icol + 3] = (
                    I6_4 * Mxyz[irow, jcol + 4]
                    + I6_11 * Mxyz[irow, jcol + 11]
                    + I6_13 * Mxyz[irow, jcol + 13]
                )
                Mtemp[irow, icol + 3] += (
                    I6_22 * Mxyz[irow, jcol + 22] + I6_24 * Mxyz[irow, jcol + 24]
                )
                # ml: -2
                Mtemp[irow, icol + 4] = (
                    I4_1 * Mxyz[irow, jcol + 1]
                    + I4_6 * Mxyz[irow, jcol + 6]
                    + I4_8 * Mxyz[irow, jcol + 8]
                )
                Mtemp[irow, icol + 4] += (
                    I4_15 * Mxyz[irow, jcol + 15]
                    + I4_17 * Mxyz[irow, jcol + 17]
                    + I4_19 * Mxyz[irow, jcol + 19]
                )
                # ml: -1
                Mtemp[irow, icol + 5] = (
                    I2_4 * Mxyz[irow, jcol + 4]
                    + I2_11 * Mxyz[irow, jcol + 11]
                    + I2_13 * Mxyz[irow, jcol + 13]
                )
                Mtemp[irow, icol + 5] += (
                    I2_22 * Mxyz[irow, jcol + 22]
                    + I2_24 * Mxyz[irow, jcol + 24]
                    + I2_26 * Mxyz[irow, jcol + 26]
                )
                # ml: 0
                Mtemp[irow, icol + 6] = (
                    I0_0 * Mxyz[irow, jcol]
                    + I0_3 * Mxyz[irow, jcol + 3]
                    + I0_5 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 6] += (
                    I0_10 * Mxyz[irow, jcol + 10]
                    + I0_12 * Mxyz[irow, jcol + 12]
                    + I0_14 * Mxyz[irow, jcol + 14]
                )
                Mtemp[irow, icol + 6] += (
                    I0_21 * Mxyz[irow, jcol + 21]
                    + I0_23 * Mxyz[irow, jcol + 23]
                    + I0_25 * Mxyz[irow, jcol + 25]
                )
                Mtemp[irow, icol + 6] += I0_27 * Mxyz[irow, jcol + 27]
                # ml: 1
                Mtemp[irow, icol + 7] = (
                    I1_2 * Mxyz[irow, jcol + 2]
                    + I1_7 * Mxyz[irow, jcol + 7]
                    + I1_9 * Mxyz[irow, jcol + 9]
                )
                Mtemp[irow, icol + 7] += (
                    I1_16 * Mxyz[irow, jcol + 16]
                    + I1_18 * Mxyz[irow, jcol + 18]
                    + I1_20 * Mxyz[irow, jcol + 20]
                )
                # ml: 2
                Mtemp[irow, icol + 8] = (
                    I3_0 * Mxyz[irow, jcol]
                    + I3_3 * Mxyz[irow, jcol + 3]
                    + I3_5 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 8] += (
                    I3_10 * Mxyz[irow, jcol + 10]
                    + I3_14 * Mxyz[irow, jcol + 14]
                    + I3_21 * Mxyz[irow, jcol + 21]
                )
                Mtemp[irow, icol + 8] += (
                    I3_23 * Mxyz[irow, jcol + 23] + I3_25 * Mxyz[irow, jcol + 25]
                )
                # ml: 3
                Mtemp[irow, icol + 9] = (
                    I5_2 * Mxyz[irow, jcol + 2]
                    + I5_7 * Mxyz[irow, jcol + 7]
                    + I5_9 * Mxyz[irow, jcol + 9]
                )
                Mtemp[irow, icol + 9] += (
                    I5_16 * Mxyz[irow, jcol + 16] + I5_18 * Mxyz[irow, jcol + 18]
                )
                # ml: 4
                Mtemp[irow, icol + 10] = (
                    I7_0 * Mxyz[irow, jcol]
                    + I7_3 * Mxyz[irow, jcol + 3]
                    + I7_5 * Mxyz[irow, jcol + 5]
                )
                Mtemp[irow, icol + 10] += (
                    I7_10 * Mxyz[irow, jcol + 10]
                    + I7_12 * Mxyz[irow, jcol + 12]
                    + I7_21 * Mxyz[irow, jcol + 21]
                )
                Mtemp[irow, icol + 10] += I7_23 * Mxyz[irow, jcol + 23]
                # ml: 5
                Mtemp[irow, icol + 11] = (
                    I9_2 * Mxyz[irow, jcol + 2]
                    + I9_7 * Mxyz[irow, jcol + 7]
                    + I9_16 * Mxyz[irow, jcol + 16]
                )
                # ml: 6
                Mtemp[irow, icol + 12] = (
                    I11_0 * Mxyz[irow, jcol]
                    + I11_3 * Mxyz[irow, jcol + 3]
                    + I11_10 * Mxyz[irow, jcol + 10]
                )
                Mtemp[irow, icol + 12] += I11_21 * Mxyz[irow, jcol + 21]
                icol += 13
                jcol += 28

    Mrtp: np.ndarray = np.zeros((Nsph, Nsph))
    irow = 0
    jrow = 0
    for ibra in TP_A:
        if ibra == "s":
            Mrtp[irow, :] = Mtemp[jrow, :]
            irow += 1
            jrow += 1
        elif ibra == "p":
            for icol in range(Nsph):
                Mrtp[irow, icol] = Mtemp[jrow, icol]
                Mrtp[irow + 1, icol] = Mtemp[jrow + 1, icol]
                Mrtp[irow + 2, icol] = Mtemp[jrow + 2, icol]
            irow += 3
            jrow += 3
        elif ibra == "d":
            for icol in range(Nsph):
                Mrtp[irow, icol] = Mtemp[jrow + 1, icol]
                Mrtp[irow + 1, icol] = Mtemp[jrow + 4, icol]
                Mrtp[irow + 2, icol] = (
                    -0.5 * Mtemp[jrow, icol]
                    - 0.5 * Mtemp[jrow + 3, icol]
                    + Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 3, icol] = Mtemp[jrow + 2, icol]
                Mrtp[irow + 4, icol] = (
                    R3_2 * Mtemp[jrow, icol] - R3_2 * Mtemp[jrow + 3, icol]
                )
            irow += 5
            jrow += 6
        elif ibra == "f":
            for icol in range(Nsph):
                Mrtp[irow, icol] = (
                    -R18_4 * Mtemp[jrow + 1, icol] - R10_4 * Mtemp[jrow + 6, icol]
                )
                Mrtp[irow + 1, icol] = Mtemp[jrow + 4, icol]
                Mrtp[irow + 2, icol] = (
                    R30_20 * Mtemp[jrow + 1, icol]
                    + R6_4 * Mtemp[jrow + 6, icol]
                    + R30_5 * Mtemp[jrow + 8, icol]
                )
                Mrtp[irow + 3, icol] = (
                    R5_03 * Mtemp[jrow + 2, icol]
                    + R5_03 * Mtemp[jrow + 7, icol]
                    + Mtemp[jrow + 9, icol]
                )
                Mrtp[irow + 4, icol] = (
                    R6_4 * Mtemp[jrow, icol]
                    + R30_20 * Mtemp[jrow + 3, icol]
                    + R30_5 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 5, icol] = (
                    R3_2 * Mtemp[jrow + 2, icol] - R3_2 * Mtemp[jrow + 7, icol]
                )
                Mrtp[irow + 6, icol] = (
                    R10_4 * Mtemp[jrow, icol] + R18_4 * Mtemp[jrow + 3, icol]
                )
            irow += 7
            jrow += 10
        elif ibra == "g":
            for icol in range(Nsph):
                Mrtp[irow, icol] = (
                    R5_2 * Mtemp[jrow + 1, icol] - R5_2 * Mtemp[jrow + 6, icol]
                )
                Mrtp[irow + 1, icol] = (
                    -R2_D3_4 * Mtemp[jrow + 4, icol] - R10_4 * Mtemp[jrow + 11, icol]
                )

                Mrtp[irow + 2, icol] = (
                    R35_14 * Mtemp[jrow + 1, icol]
                    + R35_14 * Mtemp[jrow + 6, icol]
                    + R7_D3_7 * Mtemp[jrow + 8, icol]
                )

                Mrtp[irow + 3, icol] = (
                    R14_D3_28 * Mtemp[jrow + 4, icol]
                    + R70_D3_28 * Mtemp[jrow + 11, icol]
                    + R70_7 * Mtemp[jrow + 13, icol]
                )

                Mrtp[irow + 4, icol] = (
                    D3_8 * Mtemp[jrow, icol]
                    + R105_D3_140 * Mtemp[jrow + 3, icol]
                    + R105_D3_35 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 4, icol] += (
                    D3_8 * Mtemp[jrow + 10, icol]
                    + R105_D3_35 * Mtemp[jrow + 12, icol]
                    + Mtemp[jrow + 14, icol]
                )

                Mrtp[irow + 5, icol] = (
                    R70_D3_28 * Mtemp[jrow + 2, icol]
                    + R14_D3_28 * Mtemp[jrow + 7, icol]
                    + R70_7 * Mtemp[jrow + 9, icol]
                )
                Mrtp[irow + 6, icol] = (
                    R5_4 * Mtemp[jrow, icol]
                    + R21_D3_14 * Mtemp[jrow + 5, icol]
                    - R5_4 * Mtemp[jrow + 10, icol]
                )
                Mrtp[irow + 6, icol] -= R21_D3_14 * Mtemp[jrow + 12, icol]
                Mrtp[irow + 7, icol] = (
                    R10_4 * Mtemp[jrow + 2, icol] + R2_D3_4 * Mtemp[jrow + 7, icol]
                )
                Mrtp[irow + 8, icol] = (
                    R35_8 * Mtemp[jrow, icol]
                    + R3_D3_4 * Mtemp[jrow + 3, icol]
                    + R35_8 * Mtemp[jrow + 10, icol]
                )
            irow += 9
            jrow += 15
        elif ibra == "h":
            for icol in range(Nsph):
                # ml: -5
                Mrtp[irow, icol] = (
                    H910 * Mtemp[jrow + 1, icol]
                    + H93 * Mtemp[jrow + 6, icol]
                    + H90 * Mtemp[jrow + 15, icol]
                )
                # ml: -4
                Mrtp[irow + 1, icol] = (
                    H39 * Mtemp[jrow + 4, icol] - H39 * Mtemp[jrow + 11, icol]
                )
                # ml: -3
                Mrtp[irow + 2, icol] = (
                    H50 * Mtemp[jrow + 1, icol]
                    - H53 * Mtemp[jrow + 6, icol]
                    - H512 * Mtemp[jrow + 8, icol]
                )
                Mrtp[irow + 2, icol] = (
                    Mrtp[irow + 2, icol]
                    - H50 * Mtemp[jrow + 15, icol]
                    - H55 * Mtemp[jrow + 17, icol]
                )
                # ml: -2
                Mrtp[irow + 3, icol] = (
                    H44 * Mtemp[jrow + 4, icol]
                    + H44 * Mtemp[jrow + 11, icol]
                    + H114 * Mtemp[jrow + 13, icol]
                )
                # ml: -1
                Mrtp[irow + 4, icol] = (
                    H110 * Mtemp[jrow + 1, icol]
                    + H13 * Mtemp[jrow + 6, icol]
                    + H112 * Mtemp[jrow + 8, icol]
                )
                Mrtp[irow + 4, icol] += (
                    H10 * Mtemp[jrow + 15, icol]
                    + H15 * Mtemp[jrow + 17, icol]
                    + H114 * Mtemp[jrow + 19, icol]
                )
                # ml: 0
                Mrtp[irow + 5, icol] = (
                    H02 * Mtemp[jrow + 2, icol]
                    + H07 * Mtemp[jrow + 7, icol]
                    + H09 * Mtemp[jrow + 9, icol]
                )
                Mrtp[irow + 5, icol] += (
                    H02 * Mtemp[jrow + 16, icol]
                    + H09 * Mtemp[jrow + 18, icol]
                    + Mtemp[jrow + 20, icol]
                )
                # ml: 1
                Mrtp[irow + 6, icol] = (
                    H10 * Mtemp[jrow, icol]
                    + H13 * Mtemp[jrow + 3, icol]
                    + H15 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 6, icol] += (
                    H110 * Mtemp[jrow + 10, icol]
                    + H112 * Mtemp[jrow + 12, icol]
                    + H114 * Mtemp[jrow + 14, icol]
                )
                # ml: 2
                Mrtp[irow + 7, icol] = (
                    H32 * Mtemp[jrow + 2, icol]
                    + H39 * Mtemp[jrow + 9, icol]
                    - H32 * Mtemp[jrow + 16, icol]
                    - H39 * Mtemp[jrow + 18, icol]
                )
                # ml: 3
                Mrtp[irow + 8, icol] = (
                    H50 * Mtemp[jrow, icol]
                    + H53 * Mtemp[jrow + 3, icol]
                    + H55 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 8, icol] += (
                    H512 * Mtemp[jrow + 12, icol] - H50 * Mtemp[jrow + 10, icol]
                )
                # ml: 4
                Mrtp[irow + 9, icol] = (
                    H72 * Mtemp[jrow + 2, icol]
                    + H77 * Mtemp[jrow + 7, icol]
                    + H72 * Mtemp[jrow + 16, icol]
                )
                # ml: 5
                Mrtp[irow + 10, icol] = (
                    H90 * Mtemp[jrow, icol]
                    + H93 * Mtemp[jrow + 3, icol]
                    + H910 * Mtemp[jrow + 10, icol]
                )  # C55
            irow += 11
            jrow += 21
        if iket == "i":
            for icol in range(Nsph):
                # ml: -6
                Mrtp[irow, icol] = (
                    I12_1 * Mtemp[jrow + 1, icol]
                    + I12_6 * Mtemp[jrow + 6, icol]
                    + I12_15 * Mtemp[jrow + 15, icol]
                )
                # ml: -5
                Mrtp[irow + 1, icol] = (
                    I10_4 * Mtemp[jrow + 4, icol]
                    + I10_11 * Mtemp[jrow + 11, icol]
                    + I10_22 * Mtemp[jrow + 22, icol]
                )
                # ml: -4
                Mrtp[irow + 2, icol] = (
                    I8_1 * Mtemp[jrow + 1, icol]
                    + I8_8 * Mtemp[jrow + 8, icol]
                    + I8_15 * Mtemp[jrow + 15, icol]
                )
                Mrtp[irow + 2, icol] += I8_17 * Mtemp[jrow + 17, icol]
                # ml: -3
                Mrtp[irow + 3, icol] = (
                    I6_4 * Mtemp[jrow + 4, icol]
                    + I6_11 * Mtemp[jrow + 11, icol]
                    + I6_13 * Mtemp[jrow + 13, icol]
                )
                Mrtp[irow + 3, icol] += (
                    I6_22 * Mtemp[jrow + 22, icol] + I6_24 * Mtemp[jrow + 24, icol]
                )
                # ml: -2
                Mrtp[irow + 4, icol] = (
                    I4_1 * Mtemp[jrow + 1, icol]
                    + I4_6 * Mtemp[jrow + 6, icol]
                    + I4_8 * Mtemp[jrow + 8, icol]
                )
                Mrtp[irow + 4, icol] += (
                    I4_15 * Mtemp[jrow + 15, icol]
                    + I4_17 * Mtemp[jrow + 17, icol]
                    + I4_19 * Mtemp[jrow + 19, icol]
                )
                # ml: -1
                Mrtp[irow + 5, icol] = (
                    I2_4 * Mtemp[jrow + 4, icol]
                    + I2_11 * Mtemp[jrow + 11, icol]
                    + I2_13 * Mtemp[jrow + 13, icol]
                )
                Mrtp[irow + 5, icol] += (
                    I2_22 * Mtemp[jrow + 22, icol]
                    + I2_24 * Mtemp[jrow + 24, icol]
                    + I2_26 * Mtemp[jrow + 26, icol]
                )
                # ml: 0
                Mrtp[irow + 6, icol] = (
                    I0_0 * Mtemp[jrow, icol]
                    + I0_3 * Mtemp[jrow + 3, icol]
                    + I0_5 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 6, icol] += (
                    I0_10 * Mtemp[jrow + 10, icol]
                    + I0_12 * Mtemp[jrow + 12, icol]
                    + I0_14 * Mtemp[jrow + 14, icol]
                )
                Mrtp[irow + 6, icol] += (
                    I0_21 * Mtemp[jrow + 21, icol]
                    + I0_23 * Mtemp[jrow + 23, icol]
                    + I0_25 * Mtemp[jrow + 25, icol]
                )
                Mrtp[irow + 6, icol] += I0_27 * Mtemp[jrow + 27, icol]
                # ml: 1
                Mrtp[irow + 7, icol] = (
                    I1_2 * Mtemp[jrow + 2, icol]
                    + I1_7 * Mtemp[jrow + 7, icol]
                    + I1_9 * Mtemp[jrow + 9, icol]
                )
                Mrtp[irow + 7, icol] += (
                    I1_16 * Mtemp[jrow + 16, icol]
                    + I1_18 * Mtemp[jrow + 18, icol]
                    + I1_20 * Mtemp[jrow + 20, icol]
                )
                # ml: 2
                Mrtp[irow + 8, icol] = (
                    I3_0 * Mtemp[jrow, icol]
                    + I3_3 * Mtemp[jrow + 3, icol]
                    + I3_5 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 8, icol] += (
                    I3_10 * Mtemp[jrow + 10, icol]
                    + I3_14 * Mtemp[jrow + 14, icol]
                    + I3_21 * Mtemp[jrow + 21, icol]
                )
                Mrtp[irow + 8, icol] += (
                    I3_23 * Mtemp[jrow + 23, icol] + I3_25 * Mtemp[jrow + 25, icol]
                )
                # ml: 3
                Mrtp[irow + 9, icol] = (
                    I5_2 * Mtemp[jrow + 2, icol]
                    + I5_7 * Mtemp[jrow + 7, icol]
                    + I5_9 * Mtemp[jrow + 9, icol]
                )
                Mrtp[irow + 9, icol] += (
                    I5_16 * Mtemp[jrow + 16, icol] + I5_18 * Mtemp[jrow + 18, icol]
                )
                # ml: 4
                Mrtp[irow + 10, icol] = (
                    I7_0 * Mtemp[jrow, icol]
                    + I7_3 * Mtemp[jrow + 3, icol]
                    + I7_5 * Mtemp[jrow + 5, icol]
                )
                Mrtp[irow + 10, icol] += (
                    I7_10 * Mtemp[jrow + 10, icol]
                    + I7_12 * Mtemp[jrow + 12, icol]
                    + I7_21 * Mtemp[jrow + 21, icol]
                )
                Mrtp[irow + 10, icol] += I7_23 * Mtemp[jrow + 23, icol]
                # ml: 5
                Mrtp[irow + 11, icol] = (
                    I9_2 * Mtemp[jrow + 2, icol]
                    + I9_7 * Mtemp[jrow + 7, icol]
                    + I9_16 * Mtemp[jrow + 16, icol]
                )
                # ml: 6
                Mrtp[irow + 12, icol] = (
                    I11_0 * Mtemp[jrow, icol]
                    + I11_3 * Mtemp[jrow + 3, icol]
                    + I11_10 * Mtemp[jrow + 10, icol]
                )
                Mrtp[irow + 12, icol] += I11_21 * Mtemp[jrow + 21, icol]
            icol += 13
            jcol += 28

    return Mrtp
