! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
! Riemann solver in the transverse direction for the acoustics equations
! on general quadrilateral grids.

! Split asdq (= A^* \Delta q, where * = + or -)
! into down-going flux difference bmasdqb (= B^- A^* \Delta q)
! and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    implicit none

    integer, intent(in) :: maxm, meqn, mwaves, mbc, mx, imp, maux
    double precision, intent(in)      :: ql(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in)      :: qr(meqn, 1-mbc:maxm+mbc)
    double precision, intent(out)   :: asdq(meqn, 1-mbc:maxm+mbc)
    double precision, intent(out) :: bmasdq(meqn, 1-mbc:maxm+mbc)
    double precision, intent(out) :: bpasdq(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in)    :: aux1(maux,    1-mbc:maxm+mbc)
    double precision, intent(in)    :: aux2(maux,    1-mbc:maxm+mbc)
    double precision, intent(in)    :: aux3(maux,    1-mbc:maxm+mbc)

    integer :: ixy, ixtran, iytran, ilenrat, i1, i, m
    double precision :: cm, cp, zm, zz, zp, xtran, ytran, asdqt, a1, a2
    double precision :: asdqp, bpasdqp, bpasdqt, bmasdqp, bmasdqt

    if (ixy.eq.1) then
        ixtran = 4
        iytran = 5
        ilenrat = 6
    else
        ixtran = 1
        iytran = 2
        ilenrat = 3
    endif

    do i=2-mbc,mx+mbc
        if (imp.eq.1) then
            i1 = i-1
        else
            i1 = i
        endif

!       pressure component of asdq:
        asdqp = asdq(1,i)
!
!       sound speed in each row of cells:
        cm = aux1(9,i1)
        cp = aux3(9,i1)
!
!       impedances:
        zm = aux1(8,i1)
        zz = aux2(8,i1)
        zp = aux3(8,i1)

!       up-going:
        xtran = aux3(ixtran,i)
        ytran = aux3(iytran,i)

        asdqt = xtran*asdq(2,i) + ytran*asdq(3,i)
        a2 = (asdqp + zz*asdqt) / (zz+zp)
        bpasdqp = a2*zp
        bpasdqt = a2

        bpasdq(1,i) = cp * bpasdqp
        bpasdq(2,i) = cp * xtran*bpasdqt
        bpasdq(3,i) = cp * ytran*bpasdqt

!       down-going:

        xtran = aux2(ixtran,i)
        ytran = aux2(iytran,i)

        asdqt = xtran*asdq(2,i) + ytran*asdq(3,i)
        a1 = (-asdqp + zz*asdqt) / (zm + zz)
        bmasdqp = -a1*zm
        bmasdqt = a1

        bmasdq(1,i) = -cm * bmasdqp
        bmasdq(2,i) = -cm * xtran*bmasdqt
        bmasdq(3,i) = -cm * ytran*bmasdqt

!       scale by ratio of lengths:
        do m=1,3
            bmasdq(m,i) = bmasdq(m,i)*aux2(ilenrat,i)
            bpasdq(m,i) = bpasdq(m,i)*aux3(ilenrat,i)
        enddo
    enddo
    return
    end subroutine rpt2
