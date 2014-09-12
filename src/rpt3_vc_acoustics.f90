! ==================================================================
subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! ==================================================================

!     # Riemann solver in the transverse direction for the acoustics equations
!     # with varying material properties

!     # auxN(i,1,:) holds impedance Z
!     # auxN(i,2,:) holds sound speed c
!     #
!     # On input,

!     #    ql,qr is the data along some one-dimensional slice, as in rpn3
!     #         This slice is
!     #             in the x-direction if ixyz=1,
!     #             in the y-direction if ixyz=2, or
!     #             in the z-direction if ixyz=3.
!     #    asdq is an array of flux differences (A^*\Dq).
!     #         asdq(i,:) is the flux difference propagating away from
!     #         the interface between cells i-1 and i.
!     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.

!     #    ixyz indicates the direction of the original Riemann solve,
!     #         called the x-like direction in the table below:

!     #               x-like direction   y-like direction   z-like direction
!     #      ixyz=1:        x                  y                  z
!     #      ixyz=2:        y                  z                  x
!     #      ixyz=3:        z                  x                  y

!     #    icoor indicates direction in which the transverse solve should
!     #         be performed.
!     #      icoor=2: split in the y-like direction.
!     #      icoor=3: split in the z-like direction.

!     #    For example,
!     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
!     #                        bmasdq = B^-A^*\Dq,
!     #                        bpasdq = B^+A^*\Dq.
!     #
!     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
!     #                        bmasdq = C^-B^*\Dq,
!     #                        bpasdq = C^+B^*\Dq.

!     #    The parameter imp is generally needed only if aux
!     #    arrays are being used, in order to access the appropriate
!     #    variable coefficients:

!     #    imp = 1 if asdq = A^- \Dq,  the left-going flux difference
!     #          2 if asdq = A^+ \Dq, the right-going flux difference

!     #    aux2(:,:,2) is a 1d slice of the aux array along the row
!     #                 where the data ql, qr lie.
!     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
!     #                 y-like direction
!     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
!     #                z-like direction


    implicit real*8(a-h,o-z)
    dimension     ql(meqn,1-mbc:maxm+mbc)
    dimension     qr(meqn,1-mbc:maxm+mbc)
    dimension   asdq(meqn,1-mbc:maxm+mbc)
    dimension bmasdq(meqn,1-mbc:maxm+mbc)
    dimension bpasdq(meqn,1-mbc:maxm+mbc)
    dimension   aux1(maux,1-mbc:maxm+mbc,3)
    dimension   aux2(maux,1-mbc:maxm+mbc,3)
    dimension   aux3(maux,1-mbc:maxm+mbc,3)


!     # set iuvw = 2,3,4, depending on which component of q represents
!     # velocity in the transverse direction in which splitting is to
!     # be performed:
    iuvw = ixyz + icoor
    if (iuvw > 4) iuvw = iuvw-3


    do 10 i=2-mbc,mx+mbc

    !        # The flux difference asdq is split into downward moving part
    !        # traveling at speed -c relative to the medium below and
    !        # an upward moving part traveling
    !        # at speed +c relative to the medium above.
    
    !        # Note that the sum of these parts does not give all of asdq
    !        # since there is also reflection at the interfaces which decreases
    !        # the flux.
    
    !        # set impendance and sound speed in each row of cells,
    !        # by selecting appropriate values from aux arrays depending on
    !        # values of icoor, imp:

    !        # imp is used to flag whether wave is going to left or right.
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
    
        if (icoor == 2) then
        !           # transverse direction is y-like direction so
        !           # auxN(:,:,2) holds data in appropriate plane and N=(1,2,3)
        !           # for row (below,at,above) the slice of q data
            zm = aux1(1,i1,2)
            zz = aux2(1,i1,2)
            zp = aux3(1,i1,2)
            cm = aux1(2,i1,2)
            c  = aux2(2,i1,2)
            cp = aux3(2,i1,2)
        else !! (icoor .eq. 3)
        !           # transverse direction is z-like direction so
        !           # aux2(:,:,N) holds data in appropriate plane and N=(1,2,3)
        !           # for row (below,at,above) the slice of q data
            zm = aux2(1,i1,1)
            zz = aux2(1,i1,2)
            zp = aux2(1,i1,3)
            cm = aux2(2,i1,1)
            c  = aux2(2,i1,2)
            cp = aux2(2,i1,3)
        endif
    

    !        # transmitted part of down-going wave:
        a1 = (-asdq(1,i) + asdq(iuvw,i)*zz) / &
        (zm + zz)

    !        # transmitted part of up-going wave:
        a2 = (asdq(1,i) + asdq(iuvw,i)*zz) / &
        (zz + zp)
    
    !        # The down-going flux difference bmasdq is the product  -c * wave
    
        bmasdq(1,i) = cm * a1*zm
        bmasdq(2,i) = 0.d0
        bmasdq(3,i) = 0.d0
        bmasdq(4,i) = 0.d0
        bmasdq(iuvw,i) = -cm * a1
    
    !        # The up-going flux difference bpasdq is the product  c * wave
    
        bpasdq(1,i) = cp * a2*zp
        bpasdq(2,i) = 0.d0
        bpasdq(3,i) = 0.d0
        bpasdq(4,i) = 0.d0
        bpasdq(iuvw,i) = cp * a2
    
    10 END DO


    return
    end subroutine rpt3


