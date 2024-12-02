! ==================================================================
subroutine rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)
! ==================================================================

!     # Double transverse Riemann solver for the acoustics equations
!     # with varying material properties.

!     # auxN(i,1,:) holds impedance Z
!     # auxN(i,2,:) holds sound speed c
!     #
!     # On input,

!     #    ql,qr is the data along some one-dimensional slice, as in rpn3
!     #         This slice is
!     #             in the x-direction if ixyz=1,
!     #             in the y-direction if ixyz=2, or
!     #             in the z-direction if ixyz=3.

!     #    bsasdq is an array of flux differences that result from a
!     #         transverse splitting (a previous call to rpt3).
!     #         This stands for B^* A^* \Dq but could represent any of
!     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
!     #         and icoor (see below).
!     #         Moreover, each * represents either + or -, as specified by
!     #         imp and impt.

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
!     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be
!     #                        split in z into
!     #                           cmbsasdq = C^-B^*A^*\Dq,
!     #                           cpbsasdq = C^+B^*A^*\Dq.
!     #
!     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
!     #                        split in x into
!     #                           cmbsasdq = A^-C^*B^*\Dq,
!     #                           cpbsasdq = A^+C^*B^*\Dq.

!     #    The parameters imp and impt are generally needed only if aux
!     #    arrays are being used, in order to access the appropriate
!     #    variable coefficients:

!     #    imp =  1 if bsasdq = B^*A^- \Dq, a left-going flux difference
!     #           2 if bsasdq = B^*A^+ \Dq, a right-going flux difference
!     #    impt = 1 if bsasdq = B^-A^* \Dq, a down-going flux difference
!     #           2 if bsasdq = B^+A^* \Dq, an up-going flux difference

!     #    aux2(:,:,2) is a 1d slice of the aux array along the row
!     #                 where the data ql, qr lie.
!     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
!     #                 y-like direction
!     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
!     #                z-like direction


    implicit none
    !Input
    integer, intent(in)  :: ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc,3), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: bsasdq
    
    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: cmbsasdq,cpbsasdq
    
    !Local
    integer :: i,i1,iuvw
    double precision :: zm,zz,zp,cm,c,cp,a1,a2


!     # set iuvw = 2,3,4, depending on which component of q represents
!     # velocity in the transverse direction in which splitting is to
!     # be performed:
    iuvw = ixyz + icoor
    if (iuvw > 4) iuvw = iuvw-3


    do 10 i=2-mbc,mx+mbc

    !        # The flux difference bsasdq is split into downward moving part
    !        # traveling at speed -c relative to the medium below and
    !        # an upward moving part traveling at speed +c relative to
    !        # the medium above, in the double-transverse direction specified
    !        # by icoor.
    
    !        # Note that the sum of these parts does not give all of bsasdq
    !        # since there is also reflection at the interfaces which decreases
    !        # the flux.
    
    
    !        # set impendance and sound speed in each row of cells,
    !        # by selecting appropriate values from aux arrays depending on
    !        # values of icoor, imp, impt:
    
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
    
        if (impt == 1) then
        !           # bsasdq propagates to "left" in transverse direction used in rpt3
        !           # so we either use auxN(:,:,1) or aux1(:,:,N) for N=1,2,3
        !           # depending on icoor:

            if (icoor == 2) then
            !           # new double-transverse direction is y-like direction
                zm = aux1(1,i1,1)
                zz = aux2(1,i1,1)
                zp = aux3(1,i1,1)
                cm = aux1(2,i1,1)
                c =  aux2(2,i1,1)
                cp = aux3(2,i1,1)
            else !! (icoor .eq. 3)
            !           # new double-transverse direction is z-like direction
                zm = aux1(1,i1,1)
                zz = aux1(1,i1,2)
                zp = aux1(1,i1,3)
                cm = aux1(2,i1,1)
                c =  aux1(2,i1,2)
                cp = aux1(2,i1,3)
            endif
        else
        !           # bsasdq propagates to "right" in transverse direction used in rpt3
        !           # so we either use auxN(:,:,3) or aux3(:,:,N) for N=1,2,3
        !           # depending on icoor:
            if (icoor == 2) then
            !           # new double-transverse direction is y-like direction
                zm = aux1(1,i1,3)
                zz = aux2(1,i1,3)
                zp = aux3(1,i1,3)
                cm = aux1(2,i1,3)
                c =  aux2(2,i1,3)
                cp = aux3(2,i1,3)
            else !! (icoor .eq. 3)
            !           # new double-transverse direction is z-like direction
                zm = aux3(1,i1,1)
                zz = aux3(1,i1,2)
                zp = aux3(1,i1,3)
                cm = aux3(2,i1,1)
                c =  aux3(2,i1,2)
                cp = aux3(2,i1,3)
            endif
        endif
    

    !        # transmitted part of down-going wave:
        a1 = (-bsasdq(1,i) + bsasdq(iuvw,i)*zz) / &
        (zm + zz)

    !        # transmitted part of up-going wave:
        a2 = (bsasdq(1,i) + bsasdq(iuvw,i)*zz) / &
        (zz + zp)
    
    !        # The down-going flux difference cmbsasdq is the product  -c * wave
    
        cmbsasdq(1,i) = cm * a1*zm
        cmbsasdq(2,i) = 0.d0
        cmbsasdq(3,i) = 0.d0
        cmbsasdq(4,i) = 0.d0
        cmbsasdq(iuvw,i) = -cm * a1
    
    !        # The up-going flux difference cpbsasdq is the product  c * wave
    
        cpbsasdq(1,i) = cp * a2*zp
        cpbsasdq(2,i) = 0.d0
        cpbsasdq(3,i) = 0.d0
        cpbsasdq(4,i) = 0.d0
        cpbsasdq(iuvw,i) = cp * a2
    
    10 END DO


    return
    end subroutine rptt3



