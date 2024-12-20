!   =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================

!     # Riemann solver in the transverse direction for the acoustics equations.

!     # auxN(1,i) holds rho
!     # auxN(2,i) holds c
!     #  N = 1 for row below
!     #      2 for this row
!     #      3 for row above

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.

    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

    !Local
    integer :: i,i1,mu,mv
    double precision :: a1,a2
    double precision :: cm,c,cp,zm,zz,zp


    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do 20 i = 2-mbc, mx+mbc

!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides

        if (imp == 1) then
!            # asdq = amdq, moving to left
            i1 = i-1
        else
!            # asdq = apdq, moving to right
            i1 = i
        endif

!        # The flux difference asdq is split into downward moving part
!        # traveling at speed -c relative to the medium below and
!        # an upward moving part traveling
!        # at speed +c relative to the medium above.

!        # Note that the sum of these parts does not give all of asdq
!        # since there is also reflection at the interfaces which decreases
!        # the flux.

!        # sound speed in each row of cells:
        cm = aux1(2,i1)
        c = aux2(2,i1)
        cp = aux3(2,i1)

!        # impedances:
        zm = aux1(1,i1)*aux1(2,i1)
        zz = aux2(1,i1)*aux2(2,i1)
        zp = aux3(1,i1)*aux3(2,i1)

    !        # transmitted part of down-going wave:
        a1 = (zz*asdq(1,i) + asdq(mv,i)) / (zz + zm)

!        # transmitted part of up-going wave:
        a2 = (asdq(1,i)*zz - asdq(mv,i)) / (zp + zz)

    !        # The down-going flux difference bmasdq is the product  -c * wave
    
        bmasdq(1,i) = -cm * a1
        bmasdq(mu,i) = 0.d0
        bmasdq(mv,i) = -cm * a1 * zm
    
    !        # The up-going flux difference bpasdq is the product  c * wave
    
        bpasdq(1,i) = cp * a2
        bpasdq(mu,i) = 0.d0
        bpasdq(mv,i) = -cp * a2 * zp
    
    20 END DO

    return
    end subroutine rpt2
