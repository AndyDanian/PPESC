!******************************************************************
!*      Purpose: This program computes the confluent              *
!*               hypergeometric function M(a,b,x) using           *
!*               subroutine CHGM                                  *
!*      Input  : a  --- Parameter                                 *
!*               b  --- Parameter ( b <> 0,-1,-2,... )            *
!*               x  --- Argument                                  *
!*      Output:  HG --- M(a,b,x)                                  *
!*      Example:                                                  *
!*                  a       b       x          M(a,b,x)           *
!*                -----------------------------------------       *
!*                 1.5     2.0    20.0     .1208527185D+09        *
!*                 4.5     2.0    20.0     .1103561117D+12        *
!*                -1.5     2.0    20.0     .1004836854D+05        *
!*                -4.5     2.0    20.0    -.3936045244D+03        *
!*                 1.5     2.0    50.0     .8231906643D+21        *
!*                 4.5     2.0    50.0     .9310512715D+25        *
!*                -1.5     2.0    50.0     .2998660728D+16        *
!*                -4.5     2.0    50.0    -.1807604439D+13        *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions jin.ece.uiuc.edu/routines/routines.html" *
!*                                                                *
!*                              F90 Release By J-P Moreau, Paris. *
!*                                     (www.jpmoreau.fr)          *
!* http://jean-pierre.moreau.pagesperso-orange.fr/f_function2.html*
!******************************************************************
SUBROUTINE CHGM(A,B,X,HG)
!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing ג(x)
!       ===================================================
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PI=3.141592653589793D0
    A0=A
    A1=A
    X0=X
    HG=0.0D0
    IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
        HG=1.0D+300
    ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
        HG=1.0D0
    ELSE IF (A.EQ.-1.0D0) THEN
        HG=1.0D0-X/B
    ELSE IF (A.EQ.B) THEN
        HG=DEXP(X)
    ELSE IF (A-B.EQ.1.0D0) THEN
        HG=(1.0D0+X/B)*DEXP(X)
    ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
        HG=(DEXP(X)-1.0D0)/X
    ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
        M=INT(-A)
        R=1.0D0
        HG=1.0D0
        DO 10 K=1,M
            R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
    ENDIF
    IF (HG.NE.0.0D0) RETURN
    IF (X.LT.0.0D0) THEN
        A=B-A
        A0=A
        X=DABS(X)
    ENDIF
    IF (A.LT.2.0D0) NL=0
    IF (A.GE.2.0D0) THEN
        NL=1
        LA=INT(A)
        A=A-LA-1.0D0
    ENDIF
    DO 30 N=0,NL
        IF (A0.GE.2.0D0) A=A+1.0D0
        IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
            HG=1.0D0
            RG=1.0D0
            DO 15 J=1,500
                RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                HG=HG+RG
                IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
        ELSE
            CALL GAMMA(A,TA)
            CALL GAMMA(B,TB)
            XG=B-A
            CALL GAMMA(XG,TBA)
            SUM1=1.0D0
            SUM2=1.0D0
            R1=1.0D0
            R2=1.0D0
            DO 20 I=1,8
                R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                SUM1=SUM1+R1
20               SUM2=SUM2+R2
            HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
            HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
            HG=HG1+HG2
        ENDIF
25         IF (N.EQ.0) Y0=HG
        IF (N.EQ.1) Y1=HG
30      CONTINUE
    IF (A0.GE.2.0D0) THEN
        DO 35 I=1,LA-1
            HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
            Y0=Y1
            Y1=HG
35            A=A+1.0D0
    ENDIF
    IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
    A=A1
    X=X0
    RETURN
END SUBROUTINE CHGM

SUBROUTINE GAMMA(X,GA)
!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GA --- ג(x)
!       ==================================================
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION G(26)
    PI=3.141592653589793D0
    IF (X.EQ.INT(X)) THEN
        IF (X.GT.0.0D0) THEN
            GA=1.0D0
            M1=X-1
            DO 10 K=2,M1
10               GA=GA*K
        ELSE
            GA=1.0D+300
        ENDIF
    ELSE
        IF (DABS(X).GT.1.0D0) THEN
            Z=DABS(X)
            M=INT(Z)
            R=1.0D0
            DO 15 K=1,M
15               R=R*(Z-K)
            Z=Z-M
        ELSE
            Z=X
        ENDIF
        DATA G/1.0D0,0.5772156649015329D0,  &
            -0.6558780715202538D0, -0.420026350340952D-1, &
            0.1665386113822915D0,-.421977345555443D-1,    &
            -.96219715278770D-2, .72189432466630D-2,      &
            -.11651675918591D-2, -.2152416741149D-3,      &
            .1280502823882D-3, -.201348547807D-4,         &
            -.12504934821D-5, .11330272320D-5,            &
            -.2056338417D-6, .61160950D-8,                &
            .50020075D-8, -.11812746D-8,                  &
            .1043427D-9, .77823D-11,                      &
            -.36968D-11, .51D-12,                         &
            -.206D-13, -.54D-14, .14D-14, .1D-15/
        GR=G(26)
        DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
        GA=1.0D0/(GR*Z)
        IF (DABS(X).GT.1.0D0) THEN
            GA=GA*R
            IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
        ENDIF
    ENDIF
    RETURN
END SUBROUTINE GAMMA
!end of file mchgm.f90

recursive function hermite_coefficient(i, j, t, r, alpha, beta) result(eij)
        integer, intent(in) :: i,j,t
        double precision, intent(in) :: r, alpha, beta

        double precision :: p, mu, eij

        !Hermite coefficient
        p = alpha + beta
        mu = alpha * beta / p
        if ((t < 0) .or. (t > (i + j))) then
            eij = 0.0
        else if ((i == j) .and. (j == t) .and. (t == 0)) then
            ! KAB
            eij = dexp(-mu * r * r)
        else if (j == 0) then
            ! decrement index i
            eij = (&
                  (1 / (2 * p)) * hermite_coefficient(i - 1, j, t - 1, r, alpha, beta)&
                - (mu * r / alpha) * hermite_coefficient(i - 1, j, t, r, alpha, beta)&
                + (t + 1) * hermite_coefficient(i - 1, j, t + 1, r, alpha, beta)&
            )
        else
            ! decrement index j
            eij = (&
                (1 / (2 * p)) * hermite_coefficient(i, j - 1, t - 1, r, alpha, beta)&
                + (mu * r / beta) * hermite_coefficient(i, j - 1, t, r, alpha, beta)&
                + (t + 1) * hermite_coefficient(i, j - 1, t + 1, r, alpha, beta)&
            )
        end if
end function hermite_coefficient

recursive function nuclear_repulsion(t, mu, nu, n, p, PKx, PKy, PKz, Rpc) result(pot)
    integer, intent(in) :: t, mu, nu, n
    double precision, intent(in) :: p, PKx, PKy, PKz, Rpc

    double precision :: pot, x, pi, HG, a, b
    x = p * Rpc * Rpc
    pot = 0.0
    pi = 3.141592653589793D0
    if ((t == mu) .and. (mu == nu) .and. (nu == 0)) then
            a = n + 0.5
            b = n + 1.5
            CALL CHGM(a, b, -x, HG)
            !write(*,*)"a, b, -x, CHGM",a,b,-x,HG
            pot = pot + (&
            (-2.0 * p)**n&
            * HG &
            / (2.0 * n + 1.0)&
            )
    else if ((t == mu) .and. (mu == 0)) then
        if (nu > 1) then
            pot = pot + (nu - 1) * nuclear_repulsion(t, mu, nu - 2, n + 1, p, PKx, PKy, PKz, Rpc)
        endif
        pot = pot + PKz * nuclear_repulsion(t, mu, nu - 1, n + 1, p, PKx, PKy, PKz, Rpc)
    else if (t == 0) then
        if (mu > 1) then
            pot = pot + (mu - 1) * nuclear_repulsion(t, mu - 2, nu, n + 1, p, PKx, PKy, PKz, Rpc)
        endif
        pot = pot + PKy * nuclear_repulsion(t, mu - 1, nu, n + 1, p, PKx, PKy, PKz, Rpc)
    else
        if (t > 1) then
            pot = pot + (t - 1) * nuclear_repulsion(t - 2, mu, nu, n + 1, p, PKx, PKy, PKz, Rpc)
        endif
        pot = pot + PKx * nuclear_repulsion(t - 1, mu, nu, n + 1, p, PKx, PKy, PKz, Rpc)
    end if
end function nuclear_repulsion

function electron_repulsion( &
    ! Alpha and Beta centers
    i, k, m, j, l, n, alpha, beta, Ax, Ay, Az, Bx, By, Bz, &
    ! Gamma and Delta centers
    u, v, w, x, y, z, gamma, delta, Cx, Cy, Cz, Dx, Dy, Dz &
    ) result(suma)
    !Recurrence to calculate the integrates the electron repulsion
    !            int phi_i phi_j 1/r phi_k phi_l dt
    !Equation 9.9.33 from Molecular Electronic-Structure Theory. T Helgaker, et al.

    integer, intent(in) :: i, k, m, j, l, n, u, v, w, x, y, z
    double precision, intent(in) :: alpha, beta, gamma, delta
    double precision, intent(in) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz

    double precision :: suma

    double precision hermite_coefficient, nuclear_repulsion

    double precision :: p, q, pq, Px, Py, Pz, Qx, Qy, Qz, pi, RPQ
    integer :: t, mu, nu, phi, tau, theta

    p = alpha + beta  ! center one gaussian, composite p's exponents (alpha, beta)
    q = gamma + delta ! center one gaussian, composite q's exponents (gamma, delta)

    Px = alpha * Ax + beta * Bx
    Px = Px / p
    Py = alpha * Ay + beta * By
    Py = Py / p
    Pz = alpha * Az + beta * Bz
    Pz = Pz / p

    Qx = gamma * Cx + delta * Dx
    Qx = Qx / q
    Qy = gamma * Cy + delta * Dy
    Qy = Qy / q
    Qz = gamma * Cz + delta * Dz
    Qz = Qz / q

    pq = p*q/(p+q)

    RPQ = dsqrt((Px - Qx)*(Px - Qx) +  (Py - Qy)*(Py - Qy) + (Pz - Qz)*(Pz - Qz))

    suma = 0.0
    pi = 3.141592653589793D0
    do 1 t = 1, i + j + 1
        do 1 mu = 1, k + l + 1
            do 1 nu = 1, m + n + 1
                do 1 phi = 1, u + x + 1
                    do 1 tau = 1, v + y + 1
                        do 1 theta = 1, w + z + 1
                            suma = suma + (&
                                hermite_coefficient(i, j, t-1, Ax - Bx, alpha, beta)&
                                * hermite_coefficient(k, l, mu-1, Ay - By, alpha, beta)&
                                * hermite_coefficient(m, n, nu-1, Az - Bz, alpha, beta)&
                                * hermite_coefficient(u, x, phi-1, Cx - Dx, gamma, delta)&
                                * hermite_coefficient(v, y, tau-1, Cy - Dy, gamma, delta)&
                                * hermite_coefficient(w, z, theta-1, Cz - Dz, gamma, delta)&
                                * (-1)**(phi+tau+theta-3)&
                                * nuclear_repulsion(&
                                    t + phi - 2,&
                                    mu + tau - 2,&
                                    nu + theta - 2,&
                                    0,&
                                    pq,&
                                    Px - Qx,&
                                    Py - Qy,&
                                    Pz - Qz,&
                                    RPQ)&
                            )

    1 continue
    suma = suma * 2.0*pi**(2.5)/(p*q*dsqrt(p + q))
end function electron_repulsion

function null_integral(a,b,c,d,i,j,k,l,m,n,r,s,t,u,v,w,coord) result (calculate)

    integer, intent(in) :: a,b,c,d,i,j,k,l,m,n,r,s,t,u,v,w
    double precision, dimension(2,3), intent(in) :: coord

    integer :: calculate
    integer :: a_b, b_d, c_d, a1, a2
    double precision :: coord_x, coord_y, coord_z

    a_b = 0
    b_d = 0
    c_d = 0
    calculate = 1
    if ((a .eq. b) .and. (b .eq. c) .and. (c .eq. d)) then
        even_ask_lx = MOD((i + l + r + u), 2)
        even_ask_ly = MOD((j + m + s + v), 2)
        even_ask_lz = MOD((k + n + t + w), 2)
        if ((even_ask_lx .ne. 0) .or. (even_ask_ly .ne. 0) .or. (even_ask_lz .ne. 0)) then
            calculate = 0
        endif
    else
        if (a .eq. b) then
            a_b = 1
            a1 = a
            a2 = b
        endif
        if (b .eq. d) then
            b_d = 1
            a1 = b
            a2 = d
        endif
        if (c .eq. d) then
            c_d = 1
            a1 = c
            a2 = d
        endif
    endif

    if ((a_b .eq. 1) .or. (b_d .eq. 1) .or. (c_d .eq. 1)) then
        coord_x = abs(coord(a1+1,1)) + abs(coord(a2+1,1))
        coord_y = abs(coord(a1+1,2)) + abs(coord(a2+1,2))
        coord_z = abs(coord(a1+1,3)) + abs(coord(a2+1,3))
        if ((MOD((i + l + r + u), 2) .ne. 0) .and. (coord_x .eq. 0.0)) then
            calculate = 0
        endif
        if ((MOD((j + m + s + v), 2) .ne. 0) .and. (coord_y .eq. 0.0)) then
            calculate = 0
        endif
        if ((MOD((k + n + t + w), 2) .ne. 0) .and. (coord_z .eq. 0.0)) then
            calculate = 0
        endif
    endif
end function null_integral

recursive function factorial(n) result(value)
    integer, intent(in) :: n

    integer :: value
    value = 0

    if (n .le. 1) then
        value = 1
    else
        value = n * factorial(n - 1)
    endif
end function factorial

function normalization(i,j,k,alpha) result(norma)
    integer, intent(in) :: i,j,k
    double precision, intent(in) :: alpha

    integer :: factorial
    double precision :: norma, pi
    pi = 3.141592653589793D0

    norma = (2.0 * alpha / pi) ** (3.0 / 4.0)
    norma = norma *&
        dsqrt((( 8.0 * alpha ) ** (i + j + k) * factorial(i) * factorial(j) * factorial(k)) &
        / (factorial(2*i) * factorial(2*j) * factorial(2*k)))

end function normalization

!f2py3 -c -m ee_fortran ee.F90 --f90flags=-O3
subroutine i2e(ee, counter, coord, mlx, mly, mlz, center, expon, n)
    integer, intent(in) :: n
    double precision, dimension(2,3), intent(in) :: coord
    integer, dimension(n), intent(in) :: mlx, mly, mlz, center
    double precision, dimension(n), intent(in) :: expon

    double precision, intent(out) :: ee(n,n,n,n)
    integer, intent(out) :: counter

    double precision :: electron_repulsion, normalization
    integer :: p, q, u, v, m, calculate

    counter = 0
    do 1 p = 1, n
        do 1 q = 1, p
            do 1 u = 1, p
                if (u < p) then
                    m = u
                else
                    m = q
                endif
                do 1 v = 1, m

                    !write(*,*)p,q,u,v
                    !Los if siguientes su ayudan a ganar tiempo
                    calculate = null_integral(center(p),center(q),center(u),center(v),&
                                    mlx(p), mly(p), mlz(p),&
                                    mlx(q), mly(q), mlz(q),&
                                    mlx(u), mly(u), mlz(u),&
                                    mlx(v), mly(v), mlz(v),&
                                    coord)

                    if (calculate .eq. 1) then
                        counter = counter + 1
                        ee(p,q,u,v) = electron_repulsion(&
                            mlx(p),&
                            mly(p),&
                            mlz(p),&
                            mlx(q),&
                            mly(q),&
                            mlz(q),&
                            expon(p),&
                            expon(q),&
                            coord(center(p)+1,1),&
                            coord(center(p)+1,2),&
                            coord(center(p)+1,3),&
                            coord(center(q)+1,1),&
                            coord(center(q)+1,2),&
                            coord(center(q)+1,3),&
                            mlx(u),&
                            mly(u),&
                            mlz(u),&
                            mlx(v),&
                            mly(v),&
                            mlz(v),&
                            expon(u),&
                            expon(v),&
                            coord(center(u)+1,1),&
                            coord(center(u)+1,2),&
                            coord(center(u)+1,3),&
                            coord(center(v)+1,1),&
                            coord(center(v)+1,2),&
                            coord(center(v)+1,3))&
                            *normalization(mlx(p), mly(p), mlz(p), expon(p))&
                            *normalization(mlx(q), mly(q), mlz(q), expon(q))&
                            *normalization(mlx(u), mly(u), mlz(u), expon(u))&
                            *normalization(mlx(v), mly(v), mlz(v), expon(v))
                            !write(*,*)" (",p,q,"|",u,v,") : ",ee
                            ee(p,q,v,u) = ee(p,q,u,v)
                            ee(q,p,u,v) = ee(p,q,u,v)
                            ee(q,p,v,u) = ee(p,q,u,v)
                            ee(u,v,p,q) = ee(p,q,u,v)
                            ee(u,v,q,p) = ee(p,q,u,v)
                            ee(v,u,p,q) = ee(p,q,u,v)
                            ee(v,u,q,p) = ee(p,q,u,v)
                    endif
    1 continue
    write(*,*) " Cantidad de I2 calculadas ",counter
end subroutine i2e