SUBROUTINE EEMF(EHF,C,intpot,intkin,inttwo,charge,cord,nelec,iprint,natoms,nprim)
!============================================================================
!
! Calcula la energía electrónica como la matriz de fock, para
! comparar con los resultados del DALTO, para así confirma que
! está bien
!
!============================================================================
    integer*4, intent(in) :: natoms, nelec, nprim, iprint
    double precision, intent(in) :: C(nprim,nprim),intpot(natoms,nprim,nprim)
    double precision, intent(in) :: charge(natoms), cord(natoms,3)
    double precision, intent(in) :: intkin(nprim,nprim),inttwo(nprim,nprim,nprim,nprim)

    double precision, dimension(nprim), intent(out) :: EHF

    double precision, dimension(nprim,nprim) :: P,Hcore,G,F,F1,Ven,CT,FOM
    double precision EE, Vnn, R(3), DXYZ
    integer I, J, K, L, ncont

    WRITE(*,*)
    CALL HEADER
    WRITE(*,*)'Calculo de la Energía Electrónica'
    CALL HEADER
    WRITE(*,*)
    WRITE(*,'(A)')'# núcleos electrones primitivas int2c : '
    WRITE(*,'(3I4,I10)')natoms,nelec,nprim,nprim*nprim*nprim*nprim
    !Matriz Densidad = sum_i^nelec/2sum_ab C_aiC*_bi
    CT = TRANSPOSE(C)
    if (iprint.eq.10.or.iprint.gt.100) then
        write(*,*)
        write(*,*)'Matriz Coeficientes'
        call pmatrix2d(nprim,C)
    endif
    DO 10 I = 1,nprim
        DO 10 J = 1,nprim
            P(I,J) = 0.0
            DO 10 K = 1, nelec/2
                P(I,J) = P(I,J) + 2.0*CT(K,I)*CT(K,J)  !C^T*C
    10  CONTINUE
    if (iprint.eq.20.or.iprint.gt.100) then
        write(*,*)
        write(*,*)'Matriz Densidad'
        call pmatrix2d(nprim,P)
    endif
    !Core Hamiltonian
    DO 15 I = 1,nprim
        DO 15 J = 1,nprim
            Ven(I,J) = 0.0
            DO K = 1, natoms
                Ven(I,J) = Ven(I,J) + intpot(K,I,J)
            ENDDO
        Hcore(I,J) = intkin(I,J) + Ven(I,J)
    15  CONTINUE
    if (iprint.eq.30.or.iprint.gt.100) then
        write(*,*)
        write(*,*)'Matriz Hcore'
        call pmatrix2d(nprim,Hcore)
    endif
    !Matriz G
    ncont = 0
    DO 20 I = 1, nprim
        DO 20 J = 1, nprim
            G(I,J) = 0.0
            DO 20 K = 1, nprim
                DO 20 L = 1, nprim
                    G(I,J) = G(I,J) + P(K,L)*(inttwo(I,J,K,L)-0.5D0*inttwo(I,L,K,J))
    20 CONTINUE
    if (iprint.eq.40.or.iprint.gt.100) then
        write(*,*)
        write(*,*)'Matriz G'
        call pmatrix2d(nprim,G)
    endif
    !Matriz Fock
    DO 30 I = 1,nprim
        DO 30 J = 1,nprim
            F(I,J) = Hcore(I,J) + G(I,J)
            F1(i,j) = F(i,j)
    30  CONTINUE
    !FOCK
    !AO TO MO
    if (iprint.gt.2) WRITE(*,*)
    if (iprint.gt.2) WRITE(*,*) '***Auto-Valores de Energía***'
    FOM = MATMUL(CT,MATMUL(F1,C))

    !CALL AOTOMO(nprim,C,F1)
    DO i=1,nprim
        EHF(i)=FOM(i,i)
        if (iprint.gt.2) WRITE(*,*)EHF(i)
    ENDDO
    !Repulsion nuclear
    Vnn = 0.0
    DO 35 I = 1, natoms-1
        J = I + 1
        DO 35 K = J, natoms
            DO 36 L = 1, 3
                R(L) = cord(J,L) - cord(I,L)
                R(L) = R(L)*R(L)
    36    CONTINUE
    DXYZ = dsqrt(R(1)+R(2)+R(3))
    Vnn = Vnn + charge(I)*charge(J)/DXYZ
    35  CONTINUE
    !Electronic energy
    EE = 0.0
    DO 40 I = 1,nprim
        DO 40 J = 1,nprim
            EE = EE + 0.5*P(I,J)*(Hcore(I,J) + F(I,J))
    40  CONTINUE
    WRITE(*,*)
    WRITE(*,*)'========'
    WRITE(*,*)'Energía Electrónica : ',EE
    WRITE(*,*)'Repulsión Nuclear :   ',Vnn
    WRITE(*,*)'***Energía Total : ***',EE + Vnn
    !HHe+ STO-2G
    ! Electronic energy -4.208703640140
    ! Nuclear repulsion:             1.366867140514
    ! Final HF energy:              -2.841836499626
    !
    ! autovalores AO :  -1.5901846681026206  -0.8340212792922369
    ! autovalores OM :  -1.63280253 -0.17248354

    !CH4 STO-2G
    !DALTON
    ! Electronic energy :  -52.040500693036
    ! Nuclear repulsion:            13.452154170126
    ! Final HF energy:             -38.588346522910
    !
    ! autovalores OM : -10.67237865    -0.89681323    -0.51033886    -0.51033805    -0.51033676
    !                    0.75492484     0.75492730     0.75492870     0.75530382 

    !H2O STO-2G
    !DALTON
    !Electronic energy : -81.932684120310
    !Nuclear repulsion:             9.193913158539
    !Final HF energy:             -72.738770961771 
    !
    ! autovalores OM :    -19.59030002    -1.24163772    -0.60114820    -0.42217281    -0.35399801
    !          0.64261220     0.78748498 
    !

END SUBROUTINE