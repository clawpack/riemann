subroutine rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)
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
    !
    !    bsasdq is an array of flux differences that result from a
    !         transverse splitting (a previous call to rpt3).  
    !         This stands for B^* A^* \Dq but could represent any of 
    !         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
    !         and icoor (see below).
    !         Moreover, each * represents either + or -, as specified by
    !         imp and impt.
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
    !        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be 
    !                        split in z into 
    !                           cmbsasdq = C^-B^*A^*\Dq,
    !                           cpbsasdq = C^+B^*A^*\Dq.
    !
    !        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
    !                        split in x into 
    !                           cmbsasdq = A^-C^*B^*\Dq,
    !                           cpbsasdq = A^+C^*B^*\Dq.
    !
    !    The parameters imp and impt are generally needed only if aux
    !    arrays are being used, in order to access the appropriate
    !    variable coefficients.
    
    implicit none

    integer, intent(in) :: ixyz, icoor, maxm, meqn, mwaves, mbc, mx, maux, imp, impt

    double precision, intent(in)  :: ql, qr, bsasdq, aux1, aux2, aux3
    double precision, intent(out) :: cmbsasdq, cpbsasdq
    dimension       ql(meqn,1-mbc:maxm+mbc)
    dimension       qr(meqn,1-mbc:maxm+mbc)
    dimension   bsasdq(meqn,1-mbc:maxm+mbc)
    dimension cmbsasdq(meqn,1-mbc:maxm+mbc)
    dimension cpbsasdq(meqn,1-mbc:maxm+mbc)
    dimension     aux1(maux,1-mbc:maxm+mbc,3)
    dimension     aux2(maux,1-mbc:maxm+mbc,3)
    dimension     aux3(maux,1-mbc:maxm+mbc,3)

    double precision :: coeff
    common /cparam/ coeff(3)

    integer :: iuvw, i
    double precision :: sc

    ! set iuvw = 1 for u, 2 for v, 3 for w component of velocity
    ! depending on transverse direction:
    iuvw = ixyz + icoor - 1
    if (iuvw.gt.3) iuvw = iuvw-3

    ! transverse wave goes up or down with speed based on data for
    ! original normal Riemann problem.  

    do i=2-mbc,mx+mbc
        sc = coeff(iuvw) * 0.5d0*(qr(1,i-1) + ql(1,i))
        cmbsasdq(1,i) = dmin1(sc, 0.d0) * bsasdq(1,i)
        cpbsasdq(1,i) = dmax1(sc, 0.d0) * bsasdq(1,i)
    end do

    return
end subroutine rptt3
