subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
    ! Riemann solver in the transverse direction for the 3d Burgers' equation
    !    q_t  +  u*(.5*q^2)_x + v*(.5*q^2)_y + w*(.5*q^2)_z = 0
    ! where u,v,w are a given scalars, stored in the vector coeff
    ! that is set in setprob.f and passed in the common block comrp.
    !
    ! On input,
    !
    !    ql,qr is the data along some one-dimensional slice, as in rpn3
    !         This slice is
    !             in the x-direction if ixyz=1,
    !             in the y-direction if ixyz=2, or
    !             in the z-direction if ixyz=3.
    !    asdq is an array of flux differences (A^*\Dq).
    !         asdq(i,:) is the flux difference propagating away from
    !         the interface between cells i-1 and i.
    !    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
    !
    !    ixyz indicates the direction of the original Riemann solve,
    !         called the x-like direction in the table below:
    !
    !               x-like direction   y-like direction   z-like direction
    !      ixyz=1:        x                  y                  z         
    !      ixyz=2:        y                  z                  x         
    !      ixyz=3:        z                  x                  y         
    !
    !    icoor indicates direction in which the transverse solve should 
    !         be performed.
    !      icoor=2: split in the y-like direction.
    !      icoor=3: split in the z-like direction.
    !
    !    For example,
    !      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
    !                        bmasdq = B^-A^*\Dq,
    !                        bpasdq = B^+A^*\Dq.
    !
    !      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into 
    !                        bmasdq = C^-B^*\Dq,
    !                        bpasdq = C^+B^*\Dq.
    !
    !    The parameter imp is generally needed only if aux
    !    arrays are being used, in order to access the appropriate
    !    variable coefficients.

    implicit none

    integer, intent(in) :: ixyz, icoor, maxm, meqn, mwaves, mbc, mx, maux, imp

    double precision, intent(in)  :: ql, qr, asdq, aux1, aux2, aux3
    double precision, intent(out) :: bmasdq, bpasdq
    dimension     ql(meqn,1-mbc:maxm+mbc)
    dimension     qr(meqn,1-mbc:maxm+mbc)
    dimension   asdq(meqn,1-mbc:maxm+mbc)
    dimension bmasdq(meqn,1-mbc:maxm+mbc)
    dimension bpasdq(meqn,1-mbc:maxm+mbc)
    dimension   aux1(maux,1-mbc:maxm+mbc,3)
    dimension   aux2(maux,1-mbc:maxm+mbc,3)
    dimension   aux3(maux,1-mbc:maxm+mbc,3)

    double precision :: coeff
    common /cparam/ coeff(3)

    integer :: i, iuvw
    double precision :: sb

    ! set iuvw = 1 for u, 2 for v, 3 for w component of velocity
    ! depending on transverse direction:
    iuvw = ixyz + icoor - 1
    if (iuvw.gt.3) iuvw = iuvw-3

    ! transverse wave goes up or down with speed based on data for
    ! original normal Riemann problem.  

    do i=2-mbc,mx+mbc
        sb = coeff(iuvw) * 0.5d0*(qr(1,i-1) + ql(1,i))
        bmasdq(1,i) = dmin1(sb, 0.d0) * asdq(1,i)
        bpasdq(1,i) = dmax1(sb, 0.d0) * asdq(1,i)
    end do

    return
end subroutine rpt3

