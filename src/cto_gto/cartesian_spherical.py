"""
Author: Mgs Andy Zapata
"""

import numpy as np #1.19.1

def Convert_xyz_rtp(Mxyz,TP_A):
    """
    Convert cartesian to spherical integrals

    Arg:
    Mxyz [array, float]: Cartesian integrals
    TP_A [array, string]: Symbol associated with main quantum number

    Return:
    Mrtp [array, float]: Spherical integrals
    """

    Nsph        = 0 #number of spherical primitives
    Np          = 0 #number of cartesian primitives

    R3_2        = 0.866025403784439  #0.5*sqrt(3)
    R5_03       = -0.670820393249937 #-3/10*sqrt(5)
    R6_4        = -0.612372435695794 #-1/4*sqrt(6)
    R30_20      = -0.273861278752583 #-1/20*sqrt(30)
    R30_5       =  1.095445115010332 #1/5*sqrt(30)
    R10_4       =  0.790569415042095 #1/4*sqrt(10)
    R18_4       = -1.060660171779821 #-3/4*sqrt(2)
    D3_8        =  0.375 #3/8
    R105_D3_140 =  0.219577516413420 #3/140*sqrt(105)
    R105_D3_35  = -0.878310065653680 #-3/35*sqrt(105)
    R70_D3_28   = -0.896421457000795 #-3/28*sqrt(70)
    R14_D3_28   = -0.400891862868637 #-3/28*sqrt(14)
    R70_7       =  1.195228609334394 #1/7*sqrt(70)
    R5_4        = -0.559016994374947 #-1/4*sqrt(5)
    R21_D3_14   =  0.981980506061966 #3/14*sqrt(21)
    R35_14      = -0.422577127364258 #-1/14*sqrt(35)
    R7_D3_7     =  1.133893419027682 #3/7*sqrt(7)
    R2_D3_4     = -1.060660171779821 #-3/4*sqrt(2)
    R35_8       =  0.739509972887452 #1/8*sqrt(35)
    R3_D3_4     = -1.299038105676658 #-3/4*sqrt(3)
    R5_2        =  1.118033988749895 #1/2*sqrt(5)

    H02 =  0.625
    H07 =  0.36596252735569994226
    H09 = -1.0910894511799619063
    H10 =  0.48412291827592711065
    H13 =  0.21128856368212914438
    H15 = -1.2677313820927748663
    H110=  0.16137430609197570355
    H112= -0.56694670951384084082
    H114=  1.2909944487358056284   
    H32 = -0.85391256382996653193
    H39 =  1.1180339887498948482
    H44 = -0.6454972243679028142
    H50 = -0.52291251658379721749
    H53 = 0.22821773229381921394
    H55 = 0.91287092917527685576
    H512= -1.2247448713915890491
    H72 = 0.73950997288745200532
    H77 = -1.2990381056766579701
    H90 = 0.7015607600201140098
    H93 = -1.5309310892394863114
    H910= 1.169267933366856683


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

    Mtemp = np.zeros((Np,Nsph),dtype=float)
    for irow in range(Np):
        icol = 0
        jcol = 0
        for iket in TP_A:
            if iket == 's':
                Mtemp[irow,icol] = Mxyz[irow,jcol]
                icol += 1
                jcol += 1
            if iket == 'p':
                Mtemp[irow,icol] = Mxyz[irow,jcol]
                icol += 1
                Mtemp[irow,icol] = Mxyz[irow,jcol+1]
                icol += 1
                Mtemp[irow,icol] = Mxyz[irow,jcol+2]
                icol += 1
                jcol += 3
            if iket == 'd':
                Mtemp[irow,icol] = R3_2*Mxyz[irow,jcol]-R3_2*Mxyz[irow,jcol+3]
                Mtemp[irow,icol+1] = Mxyz[irow,jcol+4]
                Mtemp[irow,icol+2] = -0.5*Mxyz[irow,jcol]-0.5*Mxyz[irow,jcol+3]+Mxyz[irow,jcol+5]
                Mtemp[irow,icol+3] = Mxyz[irow,jcol+2]
                Mtemp[irow,icol+4] = Mxyz[irow,jcol+1]
                icol += 5
                jcol += 6
            if iket == 'f':
                Mtemp[irow,icol]   = -R18_4*Mxyz[irow,jcol+1]-R10_4*Mxyz[irow,jcol+6]
                Mtemp[irow,icol+1] = R3_2*Mxyz[irow,jcol+2]-R3_2*Mxyz[irow,jcol+7]
                Mtemp[irow,icol+2] = R30_20*Mxyz[irow,jcol+1]+R6_4*Mxyz[irow,jcol+6]+R30_5*Mxyz[irow,jcol+8]
                Mtemp[irow,icol+3] = R5_03*Mxyz[irow,jcol+2]+R5_03*Mxyz[irow,jcol+7]+Mxyz[irow,jcol+9]
                Mtemp[irow,icol+4] = R6_4*Mxyz[irow,jcol]+R30_20*Mxyz[irow,jcol+3]+R30_5*Mxyz[irow,jcol+5]
                Mtemp[irow,icol+5] = Mxyz[irow,jcol+4]
                Mtemp[irow,icol+6] = R10_4*Mxyz[irow,jcol]+R18_4*Mxyz[irow,jcol+3]
                icol += 7
                jcol += 10
            if iket == 'g':
                Mtemp[irow,icol]   = R35_8*Mxyz[irow,jcol]+R3_D3_4*Mxyz[irow,jcol+3]+R35_8*Mxyz[irow,jcol+10]
                Mtemp[irow,icol+1] = -R2_D3_4*Mxyz[irow,jcol+4]-R10_4*Mxyz[irow,jcol+11]

                Mtemp[irow,icol+2] = R5_4*Mxyz[irow,jcol]+R21_D3_14*Mxyz[irow,jcol+5]-R5_4*Mxyz[irow,jcol+10]
                Mtemp[irow,icol+2] -= R21_D3_14*Mxyz[irow,jcol+12]

                Mtemp[irow,icol+3] = R14_D3_28*Mxyz[irow,jcol+4]+R70_D3_28*Mxyz[irow,jcol+11]+R70_7*Mxyz[irow,jcol+13]
                
                Mtemp[irow,icol+4] = D3_8*Mxyz[irow,jcol]+R105_D3_140*Mxyz[irow,jcol+3]+R105_D3_35*Mxyz[irow,jcol+5]
                Mtemp[irow,icol+4] += D3_8*Mxyz[irow,jcol+10]+R105_D3_35*Mxyz[irow,jcol+12]+Mxyz[irow,jcol+14]

                Mtemp[irow,icol+5] = R70_D3_28*Mxyz[irow,jcol+2]+R14_D3_28*Mxyz[irow,jcol+7]+R70_7*Mxyz[irow,jcol+9]				
                Mtemp[irow,icol+6] = R35_14*Mxyz[irow,jcol+1]+R35_14*Mxyz[irow,jcol+6]+R7_D3_7*Mxyz[irow,jcol+8]
                Mtemp[irow,icol+7] = R10_4*Mxyz[irow,jcol+2]+R2_D3_4*Mxyz[irow,jcol+7]
                Mtemp[irow,icol+8] = R5_2*Mxyz[irow,jcol+1]-R5_2*Mxyz[irow,jcol+6]
                icol += 9
                jcol += 15
            if iket == 'h' :
                Mtemp[irow,icol]   = H910*Mxyz[irow,jcol+1]+H93*Mxyz[irow,jcol+6]+H90*Mxyz[irow,jcol+15]

                Mtemp[irow,icol+1] = H72*Mxyz[irow,jcol+2]+H77*Mxyz[irow,jcol+7]+H72*Mxyz[irow,jcol+16]

                Mtemp[irow,icol+2] = H50*Mxyz[irow,jcol+1]-H53*Mxyz[irow,jcol+6]-H512*Mxyz[irow,jcol+8]  #Malo
                Mtemp[irow,icol+2] = Mtemp[irow,icol+2]-H50*Mxyz[irow,jcol+15]-H55*Mxyz[irow,jcol+17]

                Mtemp[irow,icol+3] = H32*Mxyz[irow,jcol+2]+H39*Mxyz[irow,jcol+9]-H32*Mxyz[irow,jcol+16]-H39*Mxyz[irow,jcol+18] #Malo
                
                Mtemp[irow,icol+4] = H110*Mxyz[irow,jcol+1]+H13*Mxyz[irow,jcol+6]+H112*Mxyz[irow,jcol+8]
                Mtemp[irow,icol+4]+= H10*Mxyz[irow,jcol+15]+H15*Mxyz[irow,jcol+17]+H114*Mxyz[irow,jcol+19]

                Mtemp[irow,icol+5] = H02*Mxyz[irow,jcol+2]+H07*Mxyz[irow,jcol+7]+H09*Mxyz[irow,jcol+9]
                Mtemp[irow,icol+5]+= H02*Mxyz[irow,jcol+16]+H09*Mxyz[irow,jcol+18]+Mxyz[irow,jcol+20]

                Mtemp[irow,icol+6] = H10*Mxyz[irow,jcol]+H13*Mxyz[irow,jcol+3]+H15*Mxyz[irow,jcol+5]
                Mtemp[irow,icol+6]+= H110*Mxyz[irow,jcol+10]+H112*Mxyz[irow,jcol+12]+H114*Mxyz[irow,jcol+14]

                Mtemp[irow,icol+7] = H44*Mxyz[irow,jcol+4]+H44*Mxyz[irow,jcol+11]+H114*Mxyz[irow,jcol+13]

                Mtemp[irow,icol+8] = H50*Mxyz[irow,jcol]+H53*Mxyz[irow,jcol+3]+H55*Mxyz[irow,jcol+5]
                Mtemp[irow,icol+8]+= H512*Mxyz[irow,jcol+12]-H50*Mxyz[irow,jcol+10]

                Mtemp[irow,icol+9] = H39*Mxyz[irow,jcol+4]-H39*Mxyz[irow,jcol+11]

                Mtemp[irow,icol+10] = H90*Mxyz[irow,jcol]+H93*Mxyz[irow,jcol+3]+H910*Mxyz[irow,jcol+10] #C55
                icol += 11
                jcol += 21

    Mrtp = np.zeros((Nsph,Nsph),dtype=float)				
    irow = 0
    jrow = 0
    for ibra in TP_A:
        if ibra == 's':
            Mrtp[irow,:] = Mtemp[jrow,:]
            irow += 1
            jrow += 1
        elif ibra == 'p':
            Mrtp[irow,:] = Mtemp[jrow,:]
            Mrtp[irow+1,:] = Mtemp[jrow+1,:]
            Mrtp[irow+2,:] = Mtemp[jrow+2,:]
            irow += 3
            jrow += 3
        elif ibra == 'd':
            for icol in range(Nsph):
                Mrtp[irow,icol] = R3_2*Mtemp[jrow,icol]-R3_2*Mtemp[jrow+3,icol]
                Mrtp[irow+1,icol] = Mtemp[jrow+4,icol]
                Mrtp[irow+2,icol] = -0.5*Mtemp[jrow,icol]-0.5*Mtemp[jrow+3,icol]+Mtemp[jrow+5,icol]
                Mrtp[irow+3,icol] = Mtemp[jrow+2,icol]
                Mrtp[irow+4,icol] = Mtemp[jrow+1,icol]
            irow += 5
            jrow += 6	
        elif ibra == 'f':
            for icol in range(Nsph):
                Mrtp[irow,icol]   = -R18_4*Mtemp[jrow+1,icol]-R10_4*Mtemp[jrow+6,icol]
                Mrtp[irow+1,icol] = R3_2*Mtemp[jrow+2,icol]-R3_2*Mtemp[jrow+7,icol]
                Mrtp[irow+2,icol] = R30_20*Mtemp[jrow+1,icol]+R6_4*Mtemp[jrow+6,icol]+R30_5*Mtemp[jrow+8,icol]
                Mrtp[irow+3,icol] = R5_03*Mtemp[jrow+2,icol]+R5_03*Mtemp[jrow+7,icol]+Mtemp[jrow+9,icol]
                Mrtp[irow+4,icol] = R6_4*Mtemp[jrow,icol]+R30_20*Mtemp[jrow+3,icol]+R30_5*Mtemp[jrow+5,icol]
                Mrtp[irow+5,icol] = Mtemp[jrow+4,icol]
                Mrtp[irow+6,icol] = R10_4*Mtemp[jrow,icol]+R18_4*Mtemp[jrow+3,icol]
            irow += 7
            jrow += 10
        elif ibra == 'g':
            for icol in range(Nsph):
                Mrtp[irow,icol]   = R35_8*Mtemp[jrow,icol]+R3_D3_4*Mtemp[jrow+3,icol]+R35_8*Mtemp[jrow+10,icol]
                Mrtp[irow+1,icol] = -R2_D3_4*Mtemp[jrow+4,icol]-R10_4*Mtemp[jrow+11,icol]

                Mrtp[irow+2,icol] = R5_4*Mtemp[jrow,icol]+R21_D3_14*Mtemp[jrow+5,icol]-R5_4*Mtemp[jrow+10,icol]
                Mrtp[irow+2,icol] -= R21_D3_14*Mtemp[jrow+12,icol]

                Mrtp[irow+3,icol] = R14_D3_28*Mtemp[jrow+4,icol]+R70_D3_28*Mtemp[jrow+11,icol]+R70_7*Mtemp[jrow+13,icol]
                
                Mrtp[irow+4,icol] = D3_8*Mtemp[jrow,icol]+R105_D3_140*Mtemp[jrow+3,icol]+R105_D3_35*Mtemp[jrow+5,icol]
                Mrtp[irow+4,icol] += D3_8*Mtemp[jrow+10,icol]+R105_D3_35*Mtemp[jrow+12,icol]+Mtemp[jrow+14,icol]

                Mrtp[irow+5,icol] = R70_D3_28*Mtemp[jrow+2,icol]+R14_D3_28*Mtemp[jrow+7,icol]+R70_7*Mtemp[jrow+9,icol]				
                Mrtp[irow+6,icol] = R35_14*Mtemp[jrow+1,icol]+R35_14*Mtemp[jrow+6,icol]+R7_D3_7*Mtemp[jrow+8,icol]
                Mrtp[irow+7,icol] = R10_4*Mtemp[jrow+2,icol]+R2_D3_4*Mtemp[jrow+7,icol]
                Mrtp[irow+8,icol] = R5_2*Mtemp[jrow+1,icol]-R5_2*Mtemp[jrow+6,icol]
            irow += 9
            jrow += 15
        elif ibra == 'h' :
            for icol in range(Nsph):
                Mrtp[irow,icol]   = H910*Mtemp[jrow+1,icol]+H93*Mtemp[jrow+6,icol]+H90*Mtemp[jrow+15,icol]

                Mrtp[irow+1,icol] = H72*Mtemp[jrow+2,icol]+H77*Mtemp[jrow+7,icol]+H72*Mtemp[jrow+16,icol]

                Mrtp[irow+2,icol] = H50*Mtemp[jrow+1,icol]-H53*Mtemp[jrow+6,icol]-H512*Mtemp[jrow+8,icol]   
                Mrtp[irow+2,icol] = Mrtp[irow+2,icol]-H50*Mtemp[jrow+15,icol]-H55*Mtemp[jrow+17,icol]

                Mrtp[irow+3,icol] = H32*Mtemp[jrow+2,icol]+H39*Mtemp[jrow+9,icol]-H32*Mtemp[jrow+16,icol]-H39*Mtemp[jrow+18,icol] 
                
                Mrtp[irow+4,icol] = H110*Mtemp[jrow+1,icol]+H13*Mtemp[jrow+6,icol]+H112*Mtemp[jrow+8,icol]
                Mrtp[irow+4,icol]+= H10*Mtemp[jrow+15,icol]+H15*Mtemp[jrow+17,icol]+H114*Mtemp[jrow+19,icol]

                Mrtp[irow+5,icol] = H02*Mtemp[jrow+2,icol]+H07*Mtemp[jrow+7,icol]+H09*Mtemp[jrow+9,icol]
                Mrtp[irow+5,icol]+= H02*Mtemp[jrow+16,icol]+H09*Mtemp[jrow+18,icol]+Mtemp[jrow+20,icol]

                Mrtp[irow+6,icol] = H10*Mtemp[jrow,icol]+H13*Mtemp[jrow+3,icol]+H15*Mtemp[jrow+5,icol]
                Mrtp[irow+6,icol]+= H110*Mtemp[jrow+10,icol]+H112*Mtemp[jrow+12,icol]+H114*Mtemp[jrow+14,icol]

                Mrtp[irow+7,icol] = H44*Mtemp[jrow+4,icol]+H44*Mtemp[jrow+11,icol]+H114*Mtemp[jrow+13,icol]

                Mrtp[irow+8,icol] = H50*Mtemp[jrow,icol]+H53*Mtemp[jrow+3,icol]+H55*Mtemp[jrow+5,icol]
                Mrtp[irow+8,icol]+= H512*Mtemp[jrow+12,icol]-H50*Mtemp[jrow+10,icol]

                Mrtp[irow+9,icol] = H39*Mtemp[jrow+4,icol]-H39*Mtemp[jrow+11,icol]

                Mrtp[irow+10,icol] = H90*Mtemp[jrow,icol]+H93*Mtemp[jrow+3,icol]+H910*Mtemp[jrow+10,icol] #C55
            irow += 11
            jrow += 21
    return Mrtp