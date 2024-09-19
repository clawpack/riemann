! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
!     # Riemann solver in the transverse direction for the scalar equation

!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq
    !Local
    integer :: i
    double precision :: stran
    double precision :: u,v
    double precision :: stranm,stranp
    common /cparam/ u,v

!     # transverse wave speeds have been computed in rpn2
!     # maux=0 and aux arrays are unused in this example.

    if (ixy == 1) then
        stran = v
    else
        stran = u
    endif

    stranm = dmin1(stran, 0.d0)
    stranp = dmax1(stran, 0.d0)
          
    do 10 i = 2-mbc, mx+mbc
        bmasdq(1,i) = stranm * asdq(1,i)
        bpasdq(1,i) = stranp * asdq(1,i)
    10 END DO

    return
    end subroutine rpt2
