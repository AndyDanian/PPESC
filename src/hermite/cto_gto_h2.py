"""
Author: Mgs Andy Zapata
"""

from libint import *

def cto_gto_h2(Mxyz,TP_A):
    """
    Convert cartesian to spherical integrals of two--bodies

    Arg:
    Mxyz [array, float]: Cartesian integrals
    TP_A [array, string]: Symbol associated with main quantum number

    Return:
    Mrtp [array, float]: Spherical integrals
    """

    Nsph: int = 0 #number of spherical primitives
    Np: int = 0 #number of cartesian primitives
    start: float = time()

    for a in TP_A:
        if a == 's':
            Nsph += 1
            Np += 1
        if a == 'p':
            Nsph += 3
            Np += 3
        if a == 'd':
            Nsph += 5
            Np += 6
        if a == 'f':
            Nsph += 7
            Np += 10
        if a == 'g':
            Nsph += 9
            Np += 15
        if a == 'h':
            Nsph += 11
            Np += 21
        if a == 'i':
            Nsph += 13
            Np += 28

    Mtemp = np.zeros((Np,Np,Np,Nsph),dtype=float)
    for i in range(Np):
        for j in range(Np):
            for a in range(Np):
                b1 = 0
                b2 = 0
                for iket in TP_A:
                    if iket == 's':
                        Mtemp[i,j,a,b1] = Mxyz[i,j,a,b2]
                        b1 += 1
                        b2 += 1
                    if iket == 'p':
                        Mtemp[i,j,a,b1] = Mxyz[i,j,a,b2]
                        Mtemp[i,j,a,b1+1] = Mxyz[i,j,a,b2+1]
                        Mtemp[i,j,a,b1+2] = Mxyz[i,j,a,b2+2]
                        b1 += 3
                        b2 += 3
                    if iket == 'd':
                        Mtemp[i,j,a,b1]   = Mxyz[i,j,a,b2+1]
                        Mtemp[i,j,a,b1+1] = Mxyz[i,j,a,b2+4]
                        Mtemp[i,j,a,b1+2] = -0.5*Mxyz[i,j,a,b2]-0.5*Mxyz[i,j,a,b2+3]+Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+3] = Mxyz[i,j,a,b2+2]
                        Mtemp[i,j,a,b1+4] = R3_2*Mxyz[i,j,a,b2]-R3_2*Mxyz[i,j,a,b2+3]
                        b1 += 5
                        b2 += 6
                    if iket == 'f':
                        Mtemp[i,j,a,b1]   = -R18_4*Mxyz[i,j,a,b2+1]-R10_4*Mxyz[i,j,a,b2+6]
                        Mtemp[i,j,a,b1+1] = Mxyz[i,j,a,b2+4]
                        Mtemp[i,j,a,b1+2] = R30_20*Mxyz[i,j,a,b2+1]+R6_4*Mxyz[i,j,a,b2+6]+R30_5*Mxyz[i,j,a,b2+8]
                        Mtemp[i,j,a,b1+3] = R5_03*Mxyz[i,j,a,b2+2]+R5_03*Mxyz[i,j,a,b2+7]+Mxyz[i,j,a,b2+9]
                        Mtemp[i,j,a,b1+4] = R6_4*Mxyz[i,j,a,b2]+R30_20*Mxyz[i,j,a,b2+3]+R30_5*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+5] = R3_2*Mxyz[i,j,a,b2+2]-R3_2*Mxyz[i,j,a,b2+7]
                        Mtemp[i,j,a,b1+6] = R10_4*Mxyz[i,j,a,b2]+R18_4*Mxyz[i,j,a,b2+3]
                        b1 += 7
                        b2 += 10
                    if iket == 'g':
                        Mtemp[i,j,a,b1]   = R5_2*Mxyz[i,j,a,b2+1]-R5_2*Mxyz[i,j,a,b2+6]
                        Mtemp[i,j,a,b1+1] = -R2_D3_4*Mxyz[i,j,a,b2+4]-R10_4*Mxyz[i,j,a,b2+11]

                        Mtemp[i,j,a,b1+2] = R35_14*Mxyz[i,j,a,b2+1]+R35_14*Mxyz[i,j,a,b2+6]+R7_D3_7*Mxyz[i,j,a,b2+8]

                        Mtemp[i,j,a,b1+3] = R14_D3_28*Mxyz[i,j,a,b2+4]+R70_D3_28*Mxyz[i,j,a,b2+11]+R70_7*Mxyz[i,j,a,b2+13]
                        # ml: 0
                        Mtemp[i,j,a,b1+4] = D3_8*Mxyz[i,j,a,b2]+R105_D3_140*Mxyz[i,j,a,b2+3]+R105_D3_35*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+4] += D3_8*Mxyz[i,j,a,b2+10]+R105_D3_35*Mxyz[i,j,a,b2+12]+Mxyz[i,j,a,b2+14]

                        Mtemp[i,j,a,b1+5] = R70_D3_28*Mxyz[i,j,a,b2+2]+R14_D3_28*Mxyz[i,j,a,b2+7]+R70_7*Mxyz[i,j,a,b2+9]
                        Mtemp[i,j,a,b1+6] = R5_4*Mxyz[i,j,a,b2]+R21_D3_14*Mxyz[i,j,a,b2+5]-R5_4*Mxyz[i,j,a,b2+10]
                        Mtemp[i,j,a,b1+6] -= R21_D3_14*Mxyz[i,j,a,b2+12]
                        Mtemp[i,j,a,b1+7] = R10_4*Mxyz[i,j,a,b2+2]+R2_D3_4*Mxyz[i,j,a,b2+7]
                        Mtemp[i,j,a,b1+8] = R35_8*Mxyz[i,j,a,b2]+R3_D3_4*Mxyz[i,j,a,b2+3]+R35_8*Mxyz[i,j,a,b2+10]
                        b1 += 9
                        b2 += 15
                    if iket == 'h' :
                        # ml: -5
                        Mtemp[i,j,a,b1]   = H910*Mxyz[i,j,a,b2+1]+H93*Mxyz[i,j,a,b2+6]+H90*Mxyz[i,j,a,b2+15]
                        # ml: -4
                        Mtemp[i,j,a,b1+1] = H39*Mxyz[i,j,a,b2+4]-H39*Mxyz[i,j,a,b2+11]
                        # ml: -3
                        Mtemp[i,j,a,b1+2] = H50*Mxyz[i,j,a,b2+1]-H53*Mxyz[i,j,a,b2+6]-H512*Mxyz[i,j,a,b2+8]
                        Mtemp[i,j,a,b1+2] = Mtemp[i,j,a,b1+2]-H50*Mxyz[i,j,a,b2+15]-H55*Mxyz[i,j,a,b2+17]
                        # ml: -2
                        Mtemp[i,j,a,b1+3] = H44*Mxyz[i,j,a,b2+4]+H44*Mxyz[i,j,a,b2+11]+H114*Mxyz[i,j,a,b2+13]
                        # ml: -1
                        Mtemp[i,j,a,b1+4] = H110*Mxyz[i,j,a,b2+1]+H13*Mxyz[i,j,a,b2+6]+H112*Mxyz[i,j,a,b2+8]
                        Mtemp[i,j,a,b1+4]+= H10*Mxyz[i,j,a,b2+15]+H15*Mxyz[i,j,a,b2+17]+H114*Mxyz[i,j,a,b2+19]
                        # ml: 0
                        Mtemp[i,j,a,b1+5] = H02*Mxyz[i,j,a,b2+2]+H07*Mxyz[i,j,a,b2+7]+H09*Mxyz[i,j,a,b2+9]
                        Mtemp[i,j,a,b1+5]+= H02*Mxyz[i,j,a,b2+16]+H09*Mxyz[i,j,a,b2+18]+Mxyz[i,j,a,b2+20]
                        # ml: 1
                        Mtemp[i,j,a,b1+6] = H10*Mxyz[i,j,a,b2]+H13*Mxyz[i,j,a,b2+3]+H15*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+6]+= H110*Mxyz[i,j,a,b2+10]+H112*Mxyz[i,j,a,b2+12]+H114*Mxyz[i,j,a,b2+14]
                        # ml: 2
                        Mtemp[i,j,a,b1+7] = H32*Mxyz[i,j,a,b2+2]+H39*Mxyz[i,j,a,b2+9]-H32*Mxyz[i,j,a,b2+16]-H39*Mxyz[i,j,a,b2+18]
                        # ml: 3
                        Mtemp[i,j,a,b1+8] = H50*Mxyz[i,j,a,b2]+H53*Mxyz[i,j,a,b2+3]+H55*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+8]+= H512*Mxyz[i,j,a,b2+12]-H50*Mxyz[i,j,a,b2+10]
                        # ml: 4
                        Mtemp[i,j,a,b1+9] = H72*Mxyz[i,j,a,b2+2]+H77*Mxyz[i,j,a,b2+7]+H72*Mxyz[i,j,a,b2+16]
                        # ml: 5
                        Mtemp[i,j,a,b1+10] = H90*Mxyz[i,j,a,b2]+H93*Mxyz[i,j,a,b2+3]+H910*Mxyz[i,j,a,b2+10] #C55
                        b1 += 11
                        b2 += 21
                    if iket == 'i':
                        # ml: -6
                        Mtemp[i,j,a,b1] = I12_1*Mxyz[i,j,a,b2+1]+I12_6*Mxyz[i,j,a,b2+6]+I12_15*Mxyz[i,j,a,b2+15]
                        # ml: -5
                        Mtemp[i,j,a,b1+1] = I10_4*Mxyz[i,j,a,b2+4]+I10_11*Mxyz[i,j,a,b2+11]+I10_22*Mxyz[i,j,a,b2+22]
                        # ml: -4
                        Mtemp[i,j,a,b1+2] = I8_1*Mxyz[i,j,a,b2+1]+I8_8*Mxyz[i,j,a,b2+8]+I8_15*Mxyz[i,j,a,b2+15]
                        Mtemp[i,j,a,b1+2]+= I8_17*Mxyz[i,j,a,b2+17]
                        # ml: -3
                        Mtemp[i,j,a,b1+3] = I6_4*Mxyz[i,j,a,b2+4]+I6_11*Mxyz[i,j,a,b2+11]+I6_13*Mxyz[i,j,a,b2+13]
                        Mtemp[i,j,a,b1+3]+= I6_22*Mxyz[i,j,a,b2+22]+I6_24*Mxyz[i,j,a,b2+24]
                        # ml: -2
                        Mtemp[i,j,a,b1+4] = I4_1*Mxyz[i,j,a,b2+1]+I4_6*Mxyz[i,j,a,b2+6]+I4_8*Mxyz[i,j,a,b2+8]
                        Mtemp[i,j,a,b1+4]+= I4_15*Mxyz[i,j,a,b2+15]+I4_17*Mxyz[i,j,a,b2+17]+I4_19*Mxyz[i,j,a,b2+19]
                        # ml: -1
                        Mtemp[i,j,a,b1+5] = I2_4*Mxyz[i,j,a,b2+4]+I2_11*Mxyz[i,j,a,b2+11]+I2_13*Mxyz[i,j,a,b2+13]
                        Mtemp[i,j,a,b1+5]+= I2_22*Mxyz[i,j,a,b2+22]+I2_24*Mxyz[i,j,a,b2+24]+I2_26*Mxyz[i,j,a,b2+26]
                        # ml: 0
                        Mtemp[i,j,a,b1+6] = I0_0*Mxyz[i,j,a,b2]+I0_3*Mxyz[i,j,a,b2+3]+I0_5*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+6]+= I0_10*Mxyz[i,j,a,b2+10]+I0_12*Mxyz[i,j,a,b2+12]+I0_14*Mxyz[i,j,a,b2+14]
                        Mtemp[i,j,a,b1+6]+= I0_21*Mxyz[i,j,a,b2+21]+I0_23*Mxyz[i,j,a,b2+23]+I0_25*Mxyz[i,j,a,b2+25]
                        Mtemp[i,j,a,b1+6]+= I0_27*Mxyz[i,j,a,b2+27]
                        # ml: 1
                        Mtemp[i,j,a,b1+7] = I1_2*Mxyz[i,j,a,b2+2]+I1_7*Mxyz[i,j,a,b2+7]+I1_9*Mxyz[i,j,a,b2+9]
                        Mtemp[i,j,a,b1+7]+= I1_16*Mxyz[i,j,a,b2+16]+I1_18*Mxyz[i,j,a,b2+18]+I1_20*Mxyz[i,j,a,b2+20]
                        # ml: 2
                        Mtemp[i,j,a,b1+8] = I3_0*Mxyz[i,j,a,b2]+I3_3*Mxyz[i,j,a,b2+3]+I3_5*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+8]+= I3_10*Mxyz[i,j,a,b2+10]+I3_14*Mxyz[i,j,a,b2+14]+I3_21*Mxyz[i,j,a,b2+21]
                        Mtemp[i,j,a,b1+8]+= I3_23*Mxyz[i,j,a,b2+23]+I3_25*Mxyz[i,j,a,b2+25]
                        # ml: 3
                        Mtemp[i,j,a,b1+9] = I5_2*Mxyz[i,j,a,b2+2]+I5_7*Mxyz[i,j,a,b2+7]+I5_9*Mxyz[i,j,a,b2+9]
                        Mtemp[i,j,a,b1+9]+= I5_16*Mxyz[i,j,a,b2+16]+I5_18*Mxyz[i,j,a,b2+18]
                        # ml: 4
                        Mtemp[i,j,a,b1+10] = I7_0*Mxyz[i,j,a,b2]+I7_3*Mxyz[i,j,a,b2+3]+I7_5*Mxyz[i,j,a,b2+5]
                        Mtemp[i,j,a,b1+10]+= I7_10*Mxyz[i,j,a,b2+10]+I7_12*Mxyz[i,j,a,b2+12]+I7_21*Mxyz[i,j,a,b2+21]
                        Mtemp[i,j,a,b1+10]+= I7_23*Mxyz[i,j,a,b2+23]
                        # ml: 5
                        Mtemp[i,j,a,b1+11] = I9_2*Mxyz[i,j,a,b2+2]+I9_7*Mxyz[i,j,a,b2+7]+I9_16*Mxyz[i,j,a,b2+16]
                        # ml: 6
                        Mtemp[i,j,a,b1+12] = I11_0*Mxyz[i,j,a,b2]+I11_3*Mxyz[i,j,a,b2+3]+I11_10*Mxyz[i,j,a,b2+10]
                        Mtemp[i,j,a,b1+12]+= I11_21*Mxyz[i,j,a,b2+21]
                        b1 += 13
                        b2 += 28

    Mtemp1 = np.zeros((Np,Np,Nsph,Nsph),dtype=float)
    for i in range(Np):
        for j in range(Np):
            a1 = 0
            a2 = 0
            for ibra in TP_A:
                if ibra == 's':
                    Mtemp1[i,j,a1,:] = Mtemp[i,j,a2,:]
                    a1 += 1
                    a2 += 1
                elif ibra == 'p':
                    Mtemp1[i,j,a1,:] = Mtemp[i,j,a2,:]
                    Mtemp1[i,j,a1+1,:] = Mtemp[i,j,a2+1,:]
                    Mtemp1[i,j,a1+2,:] = Mtemp[i,j,a2+2,:]
                    a1 += 3
                    a2 += 3
                elif ibra == 'd':
                    for b in range(Nsph):
                        Mtemp1[i,j,a1,b]   = Mtemp[i,j,a2+1,b]
                        Mtemp1[i,j,a1+1,b] = Mtemp[i,j,a2+4,b]
                        Mtemp1[i,j,a1+2,b] = -0.5*Mtemp[i,j,a2,b]-0.5*Mtemp[i,j,a2+3,b]+Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+3,b] = Mtemp[i,j,a2+2,b]
                        Mtemp1[i,j,a1+4,b] = R3_2*Mtemp[i,j,a2,b]-R3_2*Mtemp[i,j,a2+3,b]
                    a1 += 5
                    a2 += 6
                elif ibra == 'f':
                    for b in range(Nsph):
                        Mtemp1[i,j,a1,b]   = -R18_4*Mtemp[i,j,a2+1,b]-R10_4*Mtemp[i,j,a2+6,b]
                        Mtemp1[i,j,a1+1,b] = Mtemp[i,j,a2+4,b]
                        Mtemp1[i,j,a1+2,b] = R30_20*Mtemp[i,j,a2+1,b]+R6_4*Mtemp[i,j,a2+6,b]+R30_5*Mtemp[i,j,a2+8,b]
                        Mtemp1[i,j,a1+3,b] = R5_03*Mtemp[i,j,a2+2,b]+R5_03*Mtemp[i,j,a2+7,b]+Mtemp[i,j,a2+9,b]
                        Mtemp1[i,j,a1+4,b] = R6_4*Mtemp[i,j,a2,b]+R30_20*Mtemp[i,j,a2+3,b]+R30_5*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+5,b] = R3_2*Mtemp[i,j,a2+2,b]-R3_2*Mtemp[i,j,a2+7,b]
                        Mtemp1[i,j,a1+6,b] = R10_4*Mtemp[i,j,a2,b]+R18_4*Mtemp[i,j,a2+3,b]
                    a1 += 7
                    a2 += 10
                elif ibra == 'g':
                    for b in range(Nsph):
                        Mtemp1[i,j,a1,b]   = R5_2*Mtemp[i,j,a2+1,b]-R5_2*Mtemp[i,j,a2+6,b]
                        Mtemp1[i,j,a1+1,b] = -R2_D3_4*Mtemp[i,j,a2+4,b]-R10_4*Mtemp[i,j,a2+11,b]

                        Mtemp1[i,j,a1+2,b] = R35_14*Mtemp[i,j,a2+1,b]+R35_14*Mtemp[i,j,a2+6,b]+R7_D3_7*Mtemp[i,j,a2+8,b]

                        Mtemp1[i,j,a1+3,b] = R14_D3_28*Mtemp[i,j,a2+4,b]+R70_D3_28*Mtemp[i,j,a2+11,b]+R70_7*Mtemp[i,j,a2+13,b]

                        Mtemp1[i,j,a1+4,b] = D3_8*Mtemp[i,j,a2,b]+R105_D3_140*Mtemp[i,j,a2+3,b]+R105_D3_35*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+4,b] += D3_8*Mtemp[i,j,a2+10,b]+R105_D3_35*Mtemp[i,j,a2+12,b]+Mtemp[i,j,a2+14,b]

                        Mtemp1[i,j,a1+5,b] = R70_D3_28*Mtemp[i,j,a2+2,b]+R14_D3_28*Mtemp[i,j,a2+7,b]+R70_7*Mtemp[i,j,a2+9,b]
                        Mtemp1[i,j,a1+6,b] = R5_4*Mtemp[i,j,a2,b]+R21_D3_14*Mtemp[i,j,a2+5,b]-R5_4*Mtemp[i,j,a2+10,b]
                        Mtemp1[i,j,a1+6,b] -= R21_D3_14*Mtemp[i,j,a2+12,b]
                        Mtemp1[i,j,a1+7,b] = R10_4*Mtemp[i,j,a2+2,b]+R2_D3_4*Mtemp[i,j,a2+7,b]
                        Mtemp1[i,j,a1+8,b] = R35_8*Mtemp[i,j,a2,b]+R3_D3_4*Mtemp[i,j,a2+3,b]+R35_8*Mtemp[i,j,a2+10,b]
                    a1 += 9
                    a2 += 15
                elif ibra == 'h' :
                    for b in range(Nsph):
                        # ml: -5
                        Mtemp1[i,j,a1,b]   = H910*Mtemp[i,j,a2+1,b]+H93*Mtemp[i,j,a2+6,b]+H90*Mtemp[i,j,a2+15,b]
                        # ml: -4
                        Mtemp1[i,j,a1+1,b] = H39*Mtemp[i,j,a2+4,b]-H39*Mtemp[i,j,a2+11,b]
                        # ml: -3
                        Mtemp1[i,j,a1+2,b] = H50*Mtemp[i,j,a2+1,b]-H53*Mtemp[i,j,a2+6,b]-H512*Mtemp[i,j,a2+8,b]
                        Mtemp1[i,j,a1+2,b] = Mtemp1[i,j,a1+2,b]-H50*Mtemp[i,j,a2+15,b]-H55*Mtemp[i,j,a2+17,b]
                        # ml: -2
                        Mtemp1[i,j,a1+3,b] = H44*Mtemp[i,j,a2+4,b]+H44*Mtemp[i,j,a2+11,b]+H114*Mtemp[i,j,a2+13,b]
                        # ml: -1
                        Mtemp1[i,j,a1+4,b] = H110*Mtemp[i,j,a2+1,b]+H13*Mtemp[i,j,a2+6,b]+H112*Mtemp[i,j,a2+8,b]
                        Mtemp1[i,j,a1+4,b]+= H10*Mtemp[i,j,a2+15,b]+H15*Mtemp[i,j,a2+17,b]+H114*Mtemp[i,j,a2+19,b]
                        # ml: 0
                        Mtemp1[i,j,a1+5,b] = H02*Mtemp[i,j,a2+2,b]+H07*Mtemp[i,j,a2+7,b]+H09*Mtemp[i,j,a2+9,b]
                        Mtemp1[i,j,a1+5,b]+= H02*Mtemp[i,j,a2+16,b]+H09*Mtemp[i,j,a2+18,b]+Mtemp[i,j,a2+20,b]
                        # ml: 1
                        Mtemp1[i,j,a1+6,b] = H10*Mtemp[i,j,a2,b]+H13*Mtemp[i,j,a2+3,b]+H15*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+6,b]+= H110*Mtemp[i,j,a2+10,b]+H112*Mtemp[i,j,a2+12,b]+H114*Mtemp[i,j,a2+14,b]
                        # ml: 2
                        Mtemp1[i,j,a1+7,b] = H32*Mtemp[i,j,a2+2,b]+H39*Mtemp[i,j,a2+9,b]-H32*Mtemp[i,j,a2+16,b]-H39*Mtemp[i,j,a2+18,b]
                        # ml: 3
                        Mtemp1[i,j,a1+8,b] = H50*Mtemp[i,j,a2,b]+H53*Mtemp[i,j,a2+3,b]+H55*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+8,b]+= H512*Mtemp[i,j,a2+12,b]-H50*Mtemp[i,j,a2+10,b]
                        # ml: 4
                        Mtemp1[i,j,a1+9,b] = H72*Mtemp[i,j,a2+2,b]+H77*Mtemp[i,j,a2+7,b]+H72*Mtemp[i,j,a2+16,b]
                        # ml: 5
                        Mtemp1[i,j,a1+10,b] = H90*Mtemp[i,j,a2,b]+H93*Mtemp[i,j,a2+3,b]+H910*Mtemp[i,j,a2+10,b] #C55
                    a1 += 11
                    a2 += 21
                if iket == 'i':
                    for b in range(Nsph):
                        # ml: -6
                        Mtemp1[i,j,a1,b] = I12_1*Mtemp[i,j,a2+1,b]+I12_6*Mtemp[i,j,a2+6,b]+I12_15*Mtemp[i,j,a2+15,b]
                        # ml: -5
                        Mtemp1[i,j,a1+1,b] = I10_4*Mtemp[i,j,a2+4,b]+I10_11*Mtemp[i,j,a2+11,b]+I10_22*Mtemp[i,j,a2+22,b]
                        # ml: -4
                        Mtemp1[i,j,a1+2,b] = I8_1*Mtemp[i,j,a2+1,b]+I8_8*Mtemp[i,j,a2+8,b]+I8_15*Mtemp[i,j,a2+15,b]
                        Mtemp1[i,j,a1+2,b]+= I8_17*Mtemp[i,j,a2+17,b]
                        # ml: -3
                        Mtemp1[i,j,a1+3,b] = I6_4*Mtemp[i,j,a2+4,b]+I6_11*Mtemp[i,j,a2+11,b]+I6_13*Mtemp[i,j,a2+13,b]
                        Mtemp1[i,j,a1+3,b]+= I6_22*Mtemp[i,j,a2+22,b]+I6_24*Mtemp[i,j,a2+24,b]
                        # ml: -2
                        Mtemp1[i,j,a1+4,b] = I4_1*Mtemp[i,j,a2+1,b]+I4_6*Mtemp[i,j,a2+6,b]+I4_8*Mtemp[i,j,a2+8,b]
                        Mtemp1[i,j,a1+4,b]+= I4_15*Mtemp[i,j,a2+15,b]+I4_17*Mtemp[i,j,a2+17,b]+I4_19*Mtemp[i,j,a2+19,b]
                        # ml: -1
                        Mtemp1[i,j,a1+5,b] = I2_4*Mtemp[i,j,a2+4,b]+I2_11*Mtemp[i,j,a2+11,b]+I2_13*Mtemp[i,j,a2+13,b]
                        Mtemp1[i,j,a1+5,b]+= I2_22*Mtemp[i,j,a2+22,b]+I2_24*Mtemp[i,j,a2+24,b]+I2_26*Mtemp[i,j,a2+26,b]
                        # ml: 0
                        Mtemp1[i,j,a1+6,b] = I0_0*Mtemp[i,j,a2,b]+I0_3*Mtemp[i,j,a2+3,b]+I0_5*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+6,b]+= I0_10*Mtemp[i,j,a2+10,b]+I0_12*Mtemp[i,j,a2+12,b]+I0_14*Mtemp[i,j,a2+14,b]
                        Mtemp1[i,j,a1+6,b]+= I0_21*Mtemp[i,j,a2+21,b]+I0_23*Mtemp[i,j,a2+23,b]+I0_25*Mtemp[i,j,a2+25,b]
                        Mtemp1[i,j,a1+6,b]+= I0_27*Mtemp[i,j,a2+27,b]
                        # ml: 1
                        Mtemp1[i,j,a1+7,b] = I1_2*Mtemp[i,j,a2+2,b]+I1_7*Mtemp[i,j,a2+7,b]+I1_9*Mtemp[i,j,a2+9,b]
                        Mtemp1[i,j,a1+7,b]+= I1_16*Mtemp[i,j,a2+16,b]+I1_18*Mtemp[i,j,a2+18,b]+I1_20*Mtemp[i,j,a2+20,b]
                        # ml: 2
                        Mtemp1[i,j,a1+8,b] = I3_0*Mtemp[i,j,a2,b]+I3_3*Mtemp[i,j,a2+3,b]+I3_5*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+8,b]+= I3_10*Mtemp[i,j,a2+10,b]+I3_14*Mtemp[i,j,a2+14,b]+I3_21*Mtemp[i,j,a2+21,b]
                        Mtemp1[i,j,a1+8,b]+= I3_23*Mtemp[i,j,a2+23,b]+I3_25*Mtemp[i,j,a2+25,b]
                        # ml: 3
                        Mtemp1[i,j,a1+9,b] = I5_2*Mtemp[i,j,a2+2,b]+I5_7*Mtemp[i,j,a2+7,b]+I5_9*Mtemp[i,j,a2+9,b]
                        Mtemp1[i,j,a1+9,b]+= I5_16*Mtemp[i,j,a2+16,b]+I5_18*Mtemp[i,j,a2+18,b]
                        # ml: 4
                        Mtemp1[i,j,a1+10,b] = I7_0*Mtemp[i,j,a2,b]+I7_3*Mtemp[i,j,a2+3,b]+I7_5*Mtemp[i,j,a2+5,b]
                        Mtemp1[i,j,a1+10,b]+= I7_10*Mtemp[i,j,a2+10,b]+I7_12*Mtemp[i,j,a2+12,b]+I7_21*Mtemp[i,j,a2+21,b]
                        Mtemp1[i,j,a1+10,b]+= I7_23*Mtemp[i,j,a2+23,b]
                        # ml: 5
                        Mtemp1[i,j,a1+11,b] = I9_2*Mtemp[i,j,a2+2,b]+I9_7*Mtemp[i,j,a2+7,b]+I9_16*Mtemp[i,j,a2+16,b]
                        # ml: 6
                        Mtemp1[i,j,a1+12,b] = I11_0*Mtemp[i,j,a2,b]+I11_3*Mtemp[i,j,a2+3,b]+I11_10*Mtemp[i,j,a2+10,b]
                        Mtemp1[i,j,a1+12,b]+= I11_21*Mtemp[i,j,a2+21,b]
                    a1 += 13
                    a2 += 28

    Mtemp = np.zeros((Np,Nsph,Nsph,Nsph),dtype=float)
    for i in range(Np):
        j1 = 0
        j2 = 0
        for ibra in TP_A:
            if ibra == 's':
                Mtemp[i,j1,:,:] = Mtemp1[i,j2,:,:]
                j1 += 1
                j2 += 1
            elif ibra == 'p':
                Mtemp[i,j1,:,:] = Mtemp1[i,j2,:,:]
                Mtemp[i,j1+1,:,:] = Mtemp1[i,j2+1,:,:]
                Mtemp[i,j1+2,:,:] = Mtemp1[i,j2+2,:,:]
                j1 += 3
                j2 += 3
            elif ibra == 'd':
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mtemp[i,j1,a,b]   = Mtemp1[i,j2+1,a,b]
                        Mtemp[i,j1+1,a,b] = Mtemp1[i,j2+4,a,b]
                        Mtemp[i,j1+2,a,b] = -0.5*Mtemp1[i,j2,a,b]-0.5*Mtemp1[i,j2+3,a,b]+Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+3,a,b] = Mtemp1[i,j2+2,a,b]
                        Mtemp[i,j1+4,a,b] = R3_2*Mtemp1[i,j2,a,b]-R3_2*Mtemp1[i,j2+3,a,b]
                j1 += 5
                j2 += 6
            elif ibra == 'f':
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mtemp[i,j1,a,b]   = -R18_4*Mtemp1[i,j2+1,a,b]-R10_4*Mtemp1[i,j2+6,a,b]
                        Mtemp[i,j1+1,a,b] = Mtemp1[i,j2+4,a,b]
                        Mtemp[i,j1+2,a,b] = R30_20*Mtemp1[i,j2+1,a,b]+R6_4*Mtemp1[i,j2+6,a,b]+R30_5*Mtemp1[i,j2+8,a,b]
                        Mtemp[i,j1+3,a,b] = R5_03*Mtemp1[i,j2+2,a,b]+R5_03*Mtemp1[i,j2+7,a,b]+Mtemp1[i,j2+9,a,b]
                        Mtemp[i,j1+4,a,b] = R6_4*Mtemp1[i,j2,a,b]+R30_20*Mtemp1[i,j2+3,a,b]+R30_5*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+5,a,b] = R3_2*Mtemp1[i,j2+2,a,b]-R3_2*Mtemp1[i,j2+7,a,b]
                        Mtemp[i,j1+6,a,b] = R10_4*Mtemp1[i,j2,a,b]+R18_4*Mtemp1[i,j2+3,a,b]
                j1 += 7
                j2 += 10
            elif ibra == 'g':
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mtemp[i,j1,a,b]   = R5_2*Mtemp1[i,j2+1,a,b]-R5_2*Mtemp1[i,j2+6,a,b]
                        Mtemp[i,j1+1,a,b] = -R2_D3_4*Mtemp1[i,j2+4,a,b]-R10_4*Mtemp1[i,j2+11,a,b]

                        Mtemp[i,j1+2,a,b] = R35_14*Mtemp1[i,j2+1,a,b]+R35_14*Mtemp1[i,j2+6,a,b]+R7_D3_7*Mtemp1[i,j2+8,a,b]

                        Mtemp[i,j1+3,a,b] = R14_D3_28*Mtemp1[i,j2+4,a,b]+R70_D3_28*Mtemp1[i,j2+11,a,b]+R70_7*Mtemp1[i,j2+13,a,b]

                        Mtemp[i,j1+4,a,b] = D3_8*Mtemp1[i,j2,a,b]+R105_D3_140*Mtemp1[i,j2+3,a,b]+R105_D3_35*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+4,a,b] += D3_8*Mtemp1[i,j2+10,a,b]+R105_D3_35*Mtemp1[i,j2+12,a,b]+Mtemp1[i,j2+14,a,b]

                        Mtemp[i,j1+5,a,b] = R70_D3_28*Mtemp1[i,j2+2,a,b]+R14_D3_28*Mtemp1[i,j2+7,a,b]+R70_7*Mtemp1[i,j2+9,a,b]
                        Mtemp[i,j1+6,a,b] = R5_4*Mtemp1[i,j2,a,b]+R21_D3_14*Mtemp1[i,j2+5,a,b]-R5_4*Mtemp1[i,j2+10,a,b]
                        Mtemp[i,j1+6,a,b] -= R21_D3_14*Mtemp1[i,j2+12,a,b]
                        Mtemp[i,j1+7,a,b] = R10_4*Mtemp1[i,j2+2,a,b]+R2_D3_4*Mtemp1[i,j2+7,a,b]
                        Mtemp[i,j1+8,a,b] = R35_8*Mtemp1[i,j2,a,b]+R3_D3_4*Mtemp1[i,j2+3,a,b]+R35_8*Mtemp1[i,j2+10,a,b]
                j1 += 9
                j2 += 15
            elif ibra == 'h' :
                for a in range(Nsph):
                    for b in range(Nsph):
                        # ml: -5
                        Mtemp[i,j1,a,b]   = H910*Mtemp1[i,j2+1,a,b]+H93*Mtemp1[i,j2+6,a,b]+H90*Mtemp1[i,j2+15,a,b]
                        # ml: -4
                        Mtemp[i,j1+1,a,b] = H39*Mtemp1[i,j2+4,a,b]-H39*Mtemp1[i,j2+11,a,b]
                        # ml: -3
                        Mtemp[i,j1+2,a,b] = H50*Mtemp1[i,j2+1,a,b]-H53*Mtemp1[i,j2+6,a,b]-H512*Mtemp1[i,j2+8,a,b]
                        Mtemp[i,j1+2,a,b] = Mtemp[i,j1+2,a,b]-H50*Mtemp1[i,j2+15,a,b]-H55*Mtemp1[i,j2+17,a,b]
                        # ml: -2
                        Mtemp[i,j1+3,a,b] = H44*Mtemp1[i,j2+4,a,b]+H44*Mtemp1[i,j2+11,a,b]+H114*Mtemp1[i,j2+13,a,b]
                        # ml: -1
                        Mtemp[i,j1+4,a,b] = H110*Mtemp1[i,j2+1,a,b]+H13*Mtemp1[i,j2+6,a,b]+H112*Mtemp1[i,j2+8,a,b]
                        Mtemp[i,j1+4,a,b]+= H10*Mtemp1[i,j2+15,a,b]+H15*Mtemp1[i,j2+17,a,b]+H114*Mtemp1[i,j2+19,a,b]
                        # ml: 0
                        Mtemp[i,j1+5,a,b] = H02*Mtemp1[i,j2+2,a,b]+H07*Mtemp1[i,j2+7,a,b]+H09*Mtemp1[i,j2+9,a,b]
                        Mtemp[i,j1+5,a,b]+= H02*Mtemp1[i,j2+16,a,b]+H09*Mtemp1[i,j2+18,a,b]+Mtemp1[i,j2+20,a,b]
                        # ml: 1
                        Mtemp[i,j1+6,a,b] = H10*Mtemp1[i,j2,a,b]+H13*Mtemp1[i,j2+3,a,b]+H15*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+6,a,b]+= H110*Mtemp1[i,j2+10,a,b]+H112*Mtemp1[i,j2+12,a,b]+H114*Mtemp1[i,j2+14,a,b]
                        # ml: 2
                        Mtemp[i,j1+7,a,b] = H32*Mtemp1[i,j2+2,a,b]+H39*Mtemp1[i,j2+9,a,b]-H32*Mtemp1[i,j2+16,a,b]-H39*Mtemp1[i,j2+18,a,b]
                        # ml: 3
                        Mtemp[i,j1+8,a,b] = H50*Mtemp1[i,j2,a,b]+H53*Mtemp1[i,j2+3,a,b]+H55*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+8,a,b]+= H512*Mtemp1[i,j2+12,a,b]-H50*Mtemp1[i,j2+10,a,b]
                        # ml: 4
                        Mtemp[i,j1+9,a,b] = H72*Mtemp1[i,j2+2,a,b]+H77*Mtemp1[i,j2+7,a,b]+H72*Mtemp1[i,j2+16,a,b]
                        # ml: 5
                        Mtemp[i,j1+10,a,b] = H90*Mtemp1[i,j2,a,b]+H93*Mtemp1[i,j2+3,a,b]+H910*Mtemp1[i,j2+10,a,b] #C55
                j1 += 11
                j2 += 21
            if iket == 'i':
                for a in range(Nsph):
                    for b in range(Nsph):
                        # ml: -6
                        Mtemp[i,j1,a,b] = I12_1*Mtemp1[i,j2+1,a,b]+I12_6*Mtemp1[i,j2+6,a,b]+I12_15*Mtemp1[i,j2+15,a,b]
                        # ml: -5
                        Mtemp[i,j1+1,a,b] = I10_4*Mtemp1[i,j2+4,a,b]+I10_11*Mtemp1[i,j2+11,a,b]+I10_22*Mtemp1[i,j2+22,a,b]
                        # ml: -4
                        Mtemp[i,j1+2,a,b] = I8_1*Mtemp1[i,j2+1,a,b]+I8_8*Mtemp1[i,j2+8,a,b]+I8_15*Mtemp1[i,j2+15,a,b]
                        Mtemp[i,j1+2,a,b]+= I8_17*Mtemp1[i,j2+17,a,b]
                        # ml: -3
                        Mtemp[i,j1+3,a,b] = I6_4*Mtemp1[i,j2+4,a,b]+I6_11*Mtemp1[i,j2+11,a,b]+I6_13*Mtemp1[i,j2+13,a,b]
                        Mtemp[i,j1+3,a,b]+= I6_22*Mtemp1[i,j2+22,a,b]+I6_24*Mtemp1[i,j2+24,a,b]
                        # ml: -2
                        Mtemp[i,j1+4,a,b] = I4_1*Mtemp1[i,j2+1,a,b]+I4_6*Mtemp1[i,j2+6,a,b]+I4_8*Mtemp1[i,j2+8,a,b]
                        Mtemp[i,j1+4,a,b]+= I4_15*Mtemp1[i,j2+15,a,b]+I4_17*Mtemp1[i,j2+17,a,b]+I4_19*Mtemp1[i,j2+19,a,b]
                        # ml: -1
                        Mtemp[i,j1+5,a,b] = I2_4*Mtemp1[i,j2+4,a,b]+I2_11*Mtemp1[i,j2+11,a,b]+I2_13*Mtemp1[i,j2+13,a,b]
                        Mtemp[i,j1+5,a,b]+= I2_22*Mtemp1[i,j2+22,a,b]+I2_24*Mtemp1[i,j2+24,a,b]+I2_26*Mtemp1[i,j2+26,a,b]
                        # ml: 0
                        Mtemp[i,j1+6,a,b] = I0_0*Mtemp1[i,j2,a,b]+I0_3*Mtemp1[i,j2+3,a,b]+I0_5*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+6,a,b]+= I0_10*Mtemp1[i,j2+10,a,b]+I0_12*Mtemp1[i,j2+12,a,b]+I0_14*Mtemp1[i,j2+14,a,b]
                        Mtemp[i,j1+6,a,b]+= I0_21*Mtemp1[i,j2+21,a,b]+I0_23*Mtemp1[i,j2+23,a,b]+I0_25*Mtemp1[i,j2+25,a,b]
                        Mtemp[i,j1+6,a,b]+= I0_27*Mtemp1[i,j2+27,a,b]
                        # ml: 1
                        Mtemp[i,j1+7,a,b] = I1_2*Mtemp1[i,j2+2,a,b]+I1_7*Mtemp1[i,j2+7,a,b]+I1_9*Mtemp1[i,j2+9,a,b]
                        Mtemp[i,j1+7,a,b]+= I1_16*Mtemp1[i,j2+16,a,b]+I1_18*Mtemp1[i,j2+18,a,b]+I1_20*Mtemp1[i,j2+20,a,b]
                        # ml: 2
                        Mtemp[i,j1+8,a,b] = I3_0*Mtemp1[i,j2,a,b]+I3_3*Mtemp1[i,j2+3,a,b]+I3_5*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+8,a,b]+= I3_10*Mtemp1[i,j2+10,a,b]+I3_14*Mtemp1[i,j2+14,a,b]+I3_21*Mtemp1[i,j2+21,a,b]
                        Mtemp[i,j1+8,a,b]+= I3_23*Mtemp1[i,j2+23,a,b]+I3_25*Mtemp1[i,j2+25,a,b]
                        # ml: 3
                        Mtemp[i,j1+9,a,b] = I5_2*Mtemp1[i,j2+2,a,b]+I5_7*Mtemp1[i,j2+7,a,b]+I5_9*Mtemp1[i,j2+9,a,b]
                        Mtemp[i,j1+9,a,b]+= I5_16*Mtemp1[i,j2+16,a,b]+I5_18*Mtemp1[i,j2+18,a,b]
                        # ml: 4
                        Mtemp[i,j1+10,a,b] = I7_0*Mtemp1[i,j2,a,b]+I7_3*Mtemp1[i,j2+3,a,b]+I7_5*Mtemp1[i,j2+5,a,b]
                        Mtemp[i,j1+10,a,b]+= I7_10*Mtemp1[i,j2+10,a,b]+I7_12*Mtemp1[i,j2+12,a,b]+I7_21*Mtemp1[i,j2+21,a,b]
                        Mtemp[i,j1+10,a,b]+= I7_23*Mtemp1[i,j2+23,a,b]
                        # ml: 5
                        Mtemp[i,j1+11,a,b] = I9_2*Mtemp1[i,j2+2,a,b]+I9_7*Mtemp1[i,j2+7,a,b]+I9_16*Mtemp1[i,j2+16,a,b]
                        # ml: 6
                        Mtemp[i,j1+12,a,b] = I11_0*Mtemp1[i,j2,a,b]+I11_3*Mtemp1[i,j2+3,a,b]+I11_10*Mtemp1[i,j2+10,a,b]
                        Mtemp[i,j1+12,a,b]+= I11_21*Mtemp1[i,j2+21,a,b]
                j1 += 13
                j2 += 28

    Mrtp = np.zeros((Nsph,Nsph,Nsph,Nsph),dtype=float)
    i1 = 0
    i2 = 0
    for ibra in TP_A:
        if ibra == 's':
            Mrtp[i1,:,:,:] = Mtemp[i2,:,:,:]
            i1 += 1
            i2 += 1
        elif ibra == 'p':
            Mrtp[i1,:,:,:] = Mtemp[i2,:,:,:]
            Mrtp[i1+1,:,:,:] = Mtemp[i2+1,:,:,:]
            Mrtp[i1+2,:,:,:] = Mtemp[i2+2,:,:,:]
            i1 += 3
            i2 += 3
        elif ibra == 'd':
            for j in range(Nsph):
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mrtp[i1,j,a,b]   = Mtemp[i2+1,j,a,b]
                        Mrtp[i1+1,j,a,b] = Mtemp[i2+4,j,a,b]
                        Mrtp[i1+2,j,a,b] = -0.5*Mtemp[i2,j,a,b]-0.5*Mtemp[i2+3,j,a,b]+Mtemp[i2+5,j,a,b]
                        Mrtp[i1+3,j,a,b] = Mtemp[i2+2,j,a,b]
                        Mrtp[i1+4,j,a,b] = R3_2*Mtemp[i2,j,a,b]-R3_2*Mtemp[i2+3,j,a,b]
            i1 += 5
            i2 += 6
        elif ibra == 'f':
            for j in range(Nsph):
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mrtp[i1,j,a,b]   = -R18_4*Mtemp[i2+1,j,a,b]-R10_4*Mtemp[i2+6,j,a,b]
                        Mrtp[i1+1,j,a,b] = Mtemp[i2+4,j,a,b]
                        Mrtp[i1+2,j,a,b] = R30_20*Mtemp[i2+1,j,a,b]+R6_4*Mtemp[i2+6,j,a,b]+R30_5*Mtemp[i2+8,j,a,b]
                        Mrtp[i1+3,j,a,b] = R5_03*Mtemp[i2+2,j,a,b]+R5_03*Mtemp[i2+7,j,a,b]+Mtemp[i2+9,j,a,b]
                        Mrtp[i1+4,j,a,b] = R6_4*Mtemp[i2,j,a,b]+R30_20*Mtemp[i2+3,j,a,b]+R30_5*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+5,j,a,b] = R3_2*Mtemp[i2+2,j,a,b]-R3_2*Mtemp[i2+7,j,a,b]
                        Mrtp[i1+6,j,a,b] = R10_4*Mtemp[i2,j,a,b]+R18_4*Mtemp[i2+3,j,a,b]
            i1 += 7
            i2 += 10
        elif ibra == 'g':
            for j in range(Nsph):
                for a in range(Nsph):
                    for b in range(Nsph):
                        Mrtp[i1,j,a,b]   = R5_2*Mtemp[i2+1,j,a,b]-R5_2*Mtemp[i2+6,j,a,b]
                        Mrtp[i1+1,j,a,b] = -R2_D3_4*Mtemp[i2+4,j,a,b]-R10_4*Mtemp[i2+11,j,a,b]

                        Mrtp[i1+2,j,a,b] = R35_14*Mtemp[i2+1,j,a,b]+R35_14*Mtemp[i2+6,j,a,b]+R7_D3_7*Mtemp[i2+8,j,a,b]

                        Mrtp[i1+3,j,a,b] = R14_D3_28*Mtemp[i2+4,j,a,b]+R70_D3_28*Mtemp[i2+11,j,a,b]+R70_7*Mtemp[i2+13,j,a,b]

                        Mrtp[i1+4,j,a,b] = D3_8*Mtemp[i2,j,a,b]+R105_D3_140*Mtemp[i2+3,j,a,b]+R105_D3_35*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+4,j,a,b] += D3_8*Mtemp[i2+10,j,a,b]+R105_D3_35*Mtemp[i2+12,j,a,b]+Mtemp[i2+14,j,a,b]

                        Mrtp[i1+5,j,a,b] = R70_D3_28*Mtemp[i2+2,j,a,b]+R14_D3_28*Mtemp[i2+7,j,a,b]+R70_7*Mtemp[i2+9,j,a,b]
                        Mrtp[i1+6,j,a,b] = R5_4*Mtemp[i2,j,a,b]+R21_D3_14*Mtemp[i2+5,j,a,b]-R5_4*Mtemp[i2+10,j,a,b]
                        Mrtp[i1+6,j,a,b] -= R21_D3_14*Mtemp[i2+12,j,a,b]
                        Mrtp[i1+7,j,a,b] = R10_4*Mtemp[i2+2,j,a,b]+R2_D3_4*Mtemp[i2+7,j,a,b]
                        Mrtp[i1+8,j,a,b] = R35_8*Mtemp[i2,j,a,b]+R3_D3_4*Mtemp[i2+3,j,a,b]+R35_8*Mtemp[i2+10,j,a,b]
            i1 += 9
            i2 += 15
        elif ibra == 'h' :
            for j in range(Nsph):
                for a in range(Nsph):
                    for b in range(Nsph):
                        # ml: -5
                        Mrtp[i1,j,a,b]   = H910*Mtemp[i2+1,j,a,b]+H93*Mtemp[i2+6,j,a,b]+H90*Mtemp[i2+15,j,a,b]
                        # ml: -4
                        Mrtp[i1+1,j,a,b] = H39*Mtemp[i2+4,j,a,b]-H39*Mtemp[i2+11,j,a,b]
                        # ml: -3
                        Mrtp[i1+2,j,a,b] = H50*Mtemp[i2+1,j,a,b]-H53*Mtemp[i2+6,j,a,b]-H512*Mtemp[i2+8,j,a,b]
                        Mrtp[i1+2,j,a,b] = Mrtp[i1+2,j,a,b]-H50*Mtemp[i2+15,j,a,b]-H55*Mtemp[i2+17,j,a,b]
                        # ml: -2
                        Mrtp[i1+3,j,a,b] = H44*Mtemp[i2+4,j,a,b]+H44*Mtemp[i2+11,j,a,b]+H114*Mtemp[i2+13,j,a,b]
                        # ml: -1
                        Mrtp[i1+4,j,a,b] = H110*Mtemp[i2+1,j,a,b]+H13*Mtemp[i2+6,j,a,b]+H112*Mtemp[i2+8,j,a,b]
                        Mrtp[i1+4,j,a,b]+= H10*Mtemp[i2+15,j,a,b]+H15*Mtemp[i2+17,j,a,b]+H114*Mtemp[i2+19,j,a,b]
                        # ml: 0
                        Mrtp[i1+5,j,a,b] = H02*Mtemp[i2+2,j,a,b]+H07*Mtemp[i2+7,j,a,b]+H09*Mtemp[i2+9,j,a,b]
                        Mrtp[i1+5,j,a,b]+= H02*Mtemp[i2+16,j,a,b]+H09*Mtemp[i2+18,j,a,b]+Mtemp[i2+20,j,a,b]
                        # ml: 1
                        Mrtp[i1+6,j,a,b] = H10*Mtemp[i2,j,a,b]+H13*Mtemp[i2+3,j,a,b]+H15*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+6,j,a,b]+= H110*Mtemp[i2+10,j,a,b]+H112*Mtemp[i2+12,j,a,b]+H114*Mtemp[i2+14,j,a,b]
                        # ml: 2
                        Mrtp[i1+7,j,a,b] = H32*Mtemp[i2+2,j,a,b]+H39*Mtemp[i2+9,j,a,b]-H32*Mtemp[i2+16,j,a,b]-H39*Mtemp[i2+18,j,a,b]
                        # ml: 3
                        Mrtp[i1+8,j,a,b] = H50*Mtemp[i2,j,a,b]+H53*Mtemp[i2+3,j,a,b]+H55*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+8,j,a,b]+= H512*Mtemp[i2+12,j,a,b]-H50*Mtemp[i2+10,j,a,b]
                        # ml: 4
                        Mrtp[i1+9,j,a,b] = H72*Mtemp[i2+2,j,a,b]+H77*Mtemp[i2+7,j,a,b]+H72*Mtemp[i2+16,j,a,b]
                        # ml: 5
                        Mrtp[i1+10,j,a,b] = H90*Mtemp[i2,j,a,b]+H93*Mtemp[i2+3,j,a,b]+H910*Mtemp[i2+10,j,a,b] #C55
            i1 += 11
            i2 += 21
        if iket == 'i':
            for j in range(Nsph):
                for a in range(Nsph):
                    for b in range(Nsph):
                        # ml: -6
                        Mrtp[i1,j,a,b] = I12_1*Mtemp[i2+1,j,a,b]+I12_6*Mtemp[i2+6,j,a,b]+I12_15*Mtemp[i2+15,j,a,b]
                        # ml: -5
                        Mrtp[i1+1,j,a,b] = I10_4*Mtemp[i2+4,j,a,b]+I10_11*Mtemp[i2+11,j,a,b]+I10_22*Mtemp[i2+22,j,a,b]
                        # ml: -4
                        Mrtp[i1+2,j,a,b] = I8_1*Mtemp[i2+1,j,a,b]+I8_8*Mtemp[i2+8,j,a,b]+I8_15*Mtemp[i2+15,j,a,b]
                        Mrtp[i1+2,j,a,b]+= I8_17*Mtemp[i2+17,j,a,b]
                        # ml: -3
                        Mrtp[i1+3,j,a,b] = I6_4*Mtemp[i2+4,j,a,b]+I6_11*Mtemp[i2+11,j,a,b]+I6_13*Mtemp[i2+13,j,a,b]
                        Mrtp[i1+3,j,a,b]+= I6_22*Mtemp[i2+22,j,a,b]+I6_24*Mtemp[i2+24,j,a,b]
                        # ml: -2
                        Mrtp[i1+4,j,a,b] = I4_1*Mtemp[i2+1,j,a,b]+I4_6*Mtemp[i2+6,j,a,b]+I4_8*Mtemp[i2+8,j,a,b]
                        Mrtp[i1+4,j,a,b]+= I4_15*Mtemp[i2+15,j,a,b]+I4_17*Mtemp[i2+17,j,a,b]+I4_19*Mtemp[i2+19,j,a,b]
                        # ml: -1
                        Mrtp[i1+5,j,a,b] = I2_4*Mtemp[i2+4,j,a,b]+I2_11*Mtemp[i2+11,j,a,b]+I2_13*Mtemp[i2+13,j,a,b]
                        Mrtp[i1+5,j,a,b]+= I2_22*Mtemp[i2+22,j,a,b]+I2_24*Mtemp[i2+24,j,a,b]+I2_26*Mtemp[i2+26,j,a,b]
                        # ml: 0
                        Mrtp[i1+6,j,a,b] = I0_0*Mtemp[i2,j,a,b]+I0_3*Mtemp[i2+3,j,a,b]+I0_5*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+6,j,a,b]+= I0_10*Mtemp[i2+10,j,a,b]+I0_12*Mtemp[i2+12,j,a,b]+I0_14*Mtemp[i2+14,j,a,b]
                        Mrtp[i1+6,j,a,b]+= I0_21*Mtemp[i2+21,j,a,b]+I0_23*Mtemp[i2+23,j,a,b]+I0_25*Mtemp[i2+25,j,a,b]
                        Mrtp[i1+6,j,a,b]+= I0_27*Mtemp[i2+27,j,a,b]
                        # ml: 1
                        Mrtp[i1+7,j,a,b] = I1_2*Mtemp[i2+2,j,a,b]+I1_7*Mtemp[i2+7,j,a,b]+I1_9*Mtemp[i2+9,j,a,b]
                        Mrtp[i1+7,j,a,b]+= I1_16*Mtemp[i2+16,j,a,b]+I1_18*Mtemp[i2+18,j,a,b]+I1_20*Mtemp[i2+20,j,a,b]
                        # ml: 2
                        Mrtp[i1+8,j,a,b] = I3_0*Mtemp[i2,j,a,b]+I3_3*Mtemp[i2+3,j,a,b]+I3_5*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+8,j,a,b]+= I3_10*Mtemp[i2+10,j,a,b]+I3_14*Mtemp[i2+14,j,a,b]+I3_21*Mtemp[i2+21,j,a,b]
                        Mrtp[i1+8,j,a,b]+= I3_23*Mtemp[i2+23,j,a,b]+I3_25*Mtemp[i2+25,j,a,b]
                        # ml: 3
                        Mrtp[i1+9,j,a,b] = I5_2*Mtemp[i2+2,j,a,b]+I5_7*Mtemp[i2+7,j,a,b]+I5_9*Mtemp[i2+9,j,a,b]
                        Mrtp[i1+9,j,a,b]+= I5_16*Mtemp[i2+16,j,a,b]+I5_18*Mtemp[i2+18,j,a,b]
                        # ml: 4
                        Mrtp[i1+10,j,a,b] = I7_0*Mtemp[i2,j,a,b]+I7_3*Mtemp[i2+3,j,a,b]+I7_5*Mtemp[i2+5,j,a,b]
                        Mrtp[i1+10,j,a,b]+= I7_10*Mtemp[i2+10,j,a,b]+I7_12*Mtemp[i2+12,j,a,b]+I7_21*Mtemp[i2+21,j,a,b]
                        Mrtp[i1+10,j,a,b]+= I7_23*Mtemp[i2+23,j,a,b]
                        # ml: 5
                        Mrtp[i1+11,j,a,b] = I9_2*Mtemp[i2+2,j,a,b]+I9_7*Mtemp[i2+7,j,a,b]+I9_16*Mtemp[i2+16,j,a,b]
                        # ml: 6
                        Mrtp[i1+12,j,a,b] = I11_0*Mtemp[i2,j,a,b]+I11_3*Mtemp[i2+3,j,a,b]+I11_10*Mtemp[i2+10,j,a,b]
                        Mrtp[i1+12,j,a,b]+= I11_21*Mtemp[i2+21,j,a,b]
            i1 += 13
            i2 += 28

    print(f"\n***Time to transform two--body cto to gto {time() - start} s\n")
    return Mrtp