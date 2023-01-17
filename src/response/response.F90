!f2py3 -c -m f90response response.F90 --f90flags="-m64 -cpp -ffixed-line-length-none -ffree-line-length-none -finit-local-zero -Ofast -mtune=native -march=native -ffast-math -mfpmath=sse -msse2 -ffast-math -g -fPIC"
subroutine ao2mo(coulomb, exchange,  &
                moco,               &
                a2i,                &
                nocc, nvir)
    integer, intent(in) :: nocc, nvir
    double precision, dimension(nocc+nvir,nocc+nvir), intent(in) :: moco
    double precision, dimension(nocc+nvir,nocc+nvir,nocc+nvir,nocc+nvir), intent(in) :: a2i

    double precision, dimension(nvir,nvir,nocc,nocc), intent(out) :: coulomb
    double precision, dimension(nvir,nocc,nvir,nocc), intent(out) :: exchange

    double precision, dimension(nocc+nvir,nocc+nvir,nocc+nvir) :: f
    double precision, dimension(nocc+nvir,nocc+nvir) :: g
    double precision, dimension(nocc+nvir) :: h
    integer :: i, j, k, l, a, s, nprim

    nprim = nocc + nvir
    do 1 a = 1, nvir
        s = a + nocc

        !Transform first index
        f = 0.0D0
        do 10 i = 1, nprim
            do 10 j = 1, nprim
                do 10 k = 1, nprim
                    do 10 l = 1, nprim
                        f(i,j,k) = f(i,j,k) + moco(s,l)*a2i(i,j,k,l)
        10 continue                

        !Coulomb
        do 20 b = a, nvir
            t = b + nocc

            ! Transform second index
            g = 0.0D0
            do 30 i = 1, nprim
                do 30 l = 1, nprim
                    do 30 k = 1, nprim
                        g(i,l) = g(i,l) + moco(t,k)*f(i,l,k)
            30 continue

            do 40 j = 1, nocc
                !Transform third index
                h = 0.0D0
                do 50 i = 1, nprim
                    do 50 k = 1, nprim
                        h(i) = h(i) + moco(j,k)*g(k,i)
                50 continue

                !Transform four-th index to get Coulomb integrals
                do 60 i = j, nocc
                    do 60 k = 1, nprim
                        coulomb(a,b,j,i) = coulomb(a,b,j,i) + moco(i,k)*h(k)
                    if (b .ge. a) then
                        coulomb(b,a,i,j) = coulomb(a,b,j,i)
                        coulomb(a,b,i,j) = coulomb(a,b,j,i)
                        coulomb(b,a,j,i) = coulomb(a,b,j,i)
                    endif
                60 continue
            40 continue
        20 continue
        !Exchange
        do 70 j = 1, nocc
            ! Transform second index
            g = 0.0D0
            do 80 i = 1, nprim
                do 80 l = 1, nprim
                    do 80 k = 1, nprim
                        g(i,l) = g(i,l) + moco(j,k)*f(i,l,k)
            80 continue

            do 90 b = a, nvir
                t  = b + nocc
                ! Transform third index
                h = 0.0D0
                do 100 i = 1, nprim
                    do 100 k = 1, nprim
                        h(i) = h(i) + moco(t,k)*g(k,i)
                100 continue

                ! Transform four-th index
                do 110 i = 1, nocc
                    do 120 k = 1, nprim
                        exchange(a,j,b,i) = exchange(a,j,b,i) + moco(i,k)*h(k)
                    120 continue
                    if (b .gt. a) exchange(b,i,a,j) = exchange(a,j,b,i)
                110 continue
            90 continue
        70 continue
    1 continue
end subroutine ao2mo

subroutine lineal_sum(xppy,       &
                      gpva,gpvb, &
                      pp,  &
                      verbose,        &
                      nocc,nvir)
    integer, intent(in) :: verbose,nocc,nvir
    double precision, dimension(2*nocc*nvir), intent(in) :: gpva, gpvb
    double precision, dimension(nocc*nvir,nocc*nvir), intent(in) :: pp

    double precision, intent(out) :: xppy

    integer :: i, ia, j, jb, a, s, b, t
    integer :: nrot, ipath, count 
    double precision :: response_value, spath

    ipath = 0
    i = 1
    a = 1
    nrot = nvir*nocc
    xppy = 0.0D0
    do 10 ia = 1, nrot
        s = a + n_mo_occ

        count = 1
        spath = 0.0E+0
        j = 1
        b = 1
        do 20 jb = 1, nrot
            t = b + n_mo_occ
            response_value = 0.5*(                   &
                                    gpva(ia)*        &
                                    pp(ia,jb)*       &
                                    gpvb(jb)         &
                                    +                &
                                    gpvb(jb+nrot)*   &
                                    pp(jb,ia)*       &
                                    gpva(ia+nrot)    &
                                 )
            xppy = xppy + response_value
            if (verbose .gt. 20 .and. count .eq. 1) then
                write(*,*)
                write(*,*)" #     i         s         t          j"
            endif
            if (verbose .gt. 20 .and. abs(response_value) .gt. 0.5) then !improve threshild to write value according its amount and basis set size
                write(*,*)count,i,s,t,j,response_value
            endif

            spath = spath + response
            count = count + 1
            ipath = ipath + 1

            b = b + 1
            if (b .gt. nvir) then
                b = 1
                j = j + 1
            endif
        20 continue
        if (verbose .gt. 20) then
            write(*,*)"----------------------------------------------------"
            write(*,*)"Total ",spath
            write(*,*)
        endif
        a = a + 1
        if (a .gt. nvir) then
            a = 1
            i = i + 1
        endif
    10 continue
end subroutine lineal_sum

subroutine quadratic_sum(xppyppz,       &
                        gpva,gpvb,gpvc, &
                        iaj,ibj,icj,    &
                        sat,sbt,sct,    &
                        ppa, ppb, ppc,  &
                        verbose,        &
                        nocc,nvir)
    integer, intent(in) :: verbose,nocc,nvir
    double precision, dimension(2*nocc*nvir), intent(in) :: gpva, gpvb, gpvc
    double precision, dimension(nocc*nocc), intent(in) :: iaj, ibj, icj
    double precision, dimension(nvir*nvir), intent(in) :: sat, sbt, sct
    double precision, dimension(nocc*nvir,nocc*nvir), intent(in) :: ppa, ppb, ppc

    double precision, intent(out) :: xppyppz

    integer :: i, ia, j, jb, a, s, b, t, c, d, cd, u, v
    integer :: nrot, ipath, count 
    double precision :: response_value, vavs_a, vavs_b, vavs_c, spath
    double precision :: appb1, appb2, appb3, appb4, appb5, appb6

    nrot = nvir*nocc

    ipath = 1
    i = 1
    a = 1
    s = 1
    do 10 ia = 1, nocc*nvir
        s = a + nocc

        count = 1
        spath = 0.0E+0
        j = 1
        b = 1
        t = 1
        do 20 jb = 1, nocc*nvir
            t = b + nocc

            u = 1
            v = 1
            c = 1
            d = 1
            do 30 cd = 1, nvir*nvir
                u = c + nocc
                v = d + nocc

                if (b .eq. c) then
                    vavs_c = icj(j+(i-1)*nocc)
                    vavs_b = ibj(j+(i-1)*nocc)
                    vavs_a = iaj(j+(i-1)*nocc)
                else
                    vavs_c = 0.0
                    vavs_b = 0.0
                    vavs_a = 0.0
                endif
                appb1 = 0.0D0; appb2 = 0.0D0; appb3 = 0.0D0
                appb4 = 0.0D0; appb5 = 0.0D0; appb6 = 0.0D0
                ! <i|A|a>P_{ia,jb}(<b|B|c>-d_{cb}<i|B|j>)P_{ic,jd}<d|C|j>
                appb1 =  gpva(ia) &
                        *ppa(ia,jb) &
                        *(sbt(c+(b-1)*nvir) - vavs_b) &
                        *ppc(c+(i-1)*nvir,d+(j-1)*nvir) &
                        *gpvc(d+(j-1)*nvir+nrot)
                !
                ! <i|B|a>P_{ia,jb}(<b|A|c>-d_{cb}<i|A|j>)P_{ic,jd}<d|C|j>
                appb3 = gpvb(ia)&
                    *ppb(ia,jb)&
                    *(sat(c+(b-1)*nvir) - vavs_a) &
                    *ppc(c+(i-1)*nvir,d+(j-1)*nvir) &
                    *gpvc(d+(j-1)*nvir+nrot)
                !
                ! <i|A|a>P_{ia,jb}(<b|C|c>-d_{cb}<i|C|j>)P_{ic,jd}<d|B|j>
                appb4 = gpva(ia)&
                    *ppa(ia,jb)&
                    *(sct(c+(b-1)*nvir) - vavs_c) &
                    *ppb(c+(i-1)*nvir,d+(j-1)*nvir) &
                    *gpvb(d+(j-1)*nvir+nrot)
                !
                ! <i|C|a>P_{ia,jb}(<b|A|c>-d_{cb}<i|A|j>)P_{ic,jd}<d|B|j>
                appb6 = gpvc(ia)&
                    *ppc(ia,jb)&
                    *(sat(c+(b-1)*nvir) - vavs_a) &
                    *ppb(c+(i-1)*nvir,d+(j-1)*nvir) &
                    *gpvb(d+(j-1)*nvir+nrot)
                !
                if (a .eq. d) then
                    vavs_c = icj(j+(i-1)*nocc)
                    vavs_b = ibj(j+(i-1)*nocc)
                else
                    vavs_c = 0.0
                    vavs_b = 0.0
                endif
                ! <i|C|c>P_{ic,jd}(<d|B|a>-d_{ad}<i|B|j>)P_{ia,jb}<b|A|j>
                appb2 = gpvc(c+(i-1)*nvir)&
                    *ppc(c+(i-1)*nvir,d+(j-1)*nvir)&
                    *(sbt(a+(d-1)*nvir) - vavs_b) &
                    *ppa(ia,jb) &
                    *gpva(jb+nrot)
                !
                ! <i|B|c>P_{ic,jd}(<d|C|a>-d_{ad}<i|C|j>)P_{ia,jb}<b|A|j>
                appb5 = gpvb(c+(i-1)*nvir)&
                    *ppb(c+(i-1)*nvir,d+(j-1)*nvir)&
                    *(sct(a+(d-1)*nvir) - vavs_c) &
                    *ppa(ia,jb) &
                    *gpva(jb+nrot)

                response_value = 0.5*(appb1 + appb2 + appb3 + appb4 + appb5 + appb6)
                xppyppz = xppyppz + response_value

                if (verbose .gt. 20 .and. count .eq. 1) then
                    write(*,*)" #     i    s    t     u    v     j"
                endif
                if (verbose .gt. 20 .and. abs(response_value) .gt. 0.1) then
                    write(*,*)count, i, s, t, u,&
                              v, j, response_value
                endif

                spath = spath + response_value
                count = count + 1
                ipath = ipath + 1                              
                d = d + 1
                if (d .gt. nvir) then
                    d = 1
                    c = c + 1
                endif
            30 continue
            b = b + 1
            if (b .gt. nvir) then
                b = 1
                j = j + 1
            endif
            20 continue
        if (verbose .gt. 20) then
            write(*,*)"----------------------------------------------------"
            Write(*,*)'Total ',spath
            write(*,*)
        endif
        a = a + 1
        if (a .gt. nvir) then
            a = 1
            i = i + 1
        endif
    10 continue
end subroutine 