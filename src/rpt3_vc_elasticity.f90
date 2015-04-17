! ==================================================================
subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! ==================================================================

!     # Riemann solver in the transverse direction for the elasticity equations
!     # with varying material properties

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


    implicit none
    integer, intent(in) :: ixyz, icoor, imp, maxm,meqn,mwaves,mbc,mx, maux
    double precision, intent(in) :: ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: asdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: aux1(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux2(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux3(maux,1-mbc:maxm+mbc,3)
    double precision, intent(out) :: bmasdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: bpasdq(meqn,1-mbc:maxm+mbc)
    integer :: i, iadj, j, sig_ii, sig_ij1, sig_ij2, sig_kk1, sig_kk2, u_i, u_j1, u_j2
    double precision :: wave(meqn,mwaves)
    double precision :: s(mwaves)
    double precision :: dsig_ii, dsig_ij1, dsig_ij2, du_i, du_j1, du_j2
    double precision :: lama, mua, bulka, cpa, csa, lamb, mub, bulkb, cpb, csb
    double precision :: lam, mu, bulk, cp, cs
    double precision :: det, a1, a2, a3, a4, a5, a6



!     # set sig_ii to the diagonal component of stress in the split direction
!     # set sig_ij1 and sig_ij2 to the off-diagonal components of stress in the split direction
!     # set u_i to the velocity component in the split direction
!     # set u_j1 and u_j2 to the velocity components orthogonal to the split direction

    if (ixyz + icoor == 5) then
        sig_ii = 1  ! sig_xx
        sig_kk1 = 2 ! sig_yy
        sig_kk2 = 3 ! sig_zz
        sig_ij1 = 4 ! sig_xy
        sig_ij2 = 5 ! sig_xz
        u_i = 7     ! u
        u_j1 = 8    ! v
        u_j2 = 9    ! w
    else if (ixyz + icoor == 3 .or. ixyz + icoor == 6) then
        sig_ii = 2  ! sig_yy
        sig_kk1 = 1 ! sig_xx
        sig_kk2 = 3 ! sig_zz
        sig_ij1 = 4 ! sig_xy
        sig_ij2 = 6 ! sig_yz
        u_i = 8     ! v
        u_j1 = 7    ! u
        u_j2 = 9    ! w
    else if (ixyz + icoor == 4) then
        sig_ii = 3  ! sig_zz
        sig_kk1 = 1 ! sig_xx
        sig_kk2 = 2 ! sig_yy
        sig_ij1 = 5 ! sig_xz
        sig_ij2 = 6 ! sig_yz
        u_i = 9     ! w
        u_j1 = 7    ! u
        u_j2 = 8    ! v
    else
        write (6,*) 'Problem in rpt3... direction mismatch'
        stop
    end if

!     # split the flux difference asdq into 3 downward parts,
!     # one traveling at speed -cp and 2 traveling at speed -cs
!     # relative to the material properties to below the interface,
!     # and 3 upward parts, one traveling at speed cp
!     # and two traveling at speed cs
!     # relative to the material properties above the interface,


    do i=2-mbc,mx+mbc

!        # imp is used to flag whether the original wave is going to left or right.
        iadj = i-2+imp    !#  =  i-1 for amdq,  i for apdq

        dsig_ii = asdq(sig_ii,i)
        dsig_ij1 = asdq(sig_ij1,i)
        dsig_ij2 = asdq(sig_ij2,i)
        du_i = asdq(u_i,i)
        du_j1 = asdq(u_j1,i)
        du_j2 = asdq(u_j2,i)

    
        if (icoor == 2) then
!           # transverse direction is y-like direction so
!           # auxN(:,:,2) holds data in appropriate plane and N=(1,2,3)
!           # for row (below,at,above) the slice of q data
            lamb = aux1(2,iadj,2)
            mub = aux1(3,iadj,2)
            bulkb = lamb + 2.d0*mub
            cpb = aux1(4,iadj,2)
            csb = aux1(5,iadj,2)

            lam = aux2(2,iadj,2)
            mu = aux2(3,iadj,2)
            bulk = lam + 2.d0*mu
            cp = aux2(4,iadj,2)
            cs = aux2(5,iadj,2)

            lama = aux3(2,iadj,2)
            mua = aux3(3,iadj,2)
            bulka = lama + 2.d0*mua
            cpa = aux3(4,iadj,2)
            csa = aux3(5,iadj,2)
        else !! (icoor .eq. 3)
!           # transverse direction is z-like direction so
!           # aux2(:,:,N) holds data in appropriate plane and N=(1,2,3)
!           # for row (below,at,above) the slice of q data
            lamb = aux2(2,iadj,1)
            mub = aux2(3,iadj,1)
            bulkb = lamb + 2.d0*mub
            cpb = aux2(4,iadj,1)
            csb = aux2(5,iadj,1)

            lam = aux2(2,iadj,2)
            mu = aux2(3,iadj,2)
            bulk = lam + 2.d0*mu
            cp = aux2(4,iadj,2)
            cs = aux2(5,iadj,2)

            lama = aux2(2,iadj,3)
            mua = aux2(3,iadj,3)
            bulka = lama + 2.d0*mua
            cpa = aux2(4,iadj,3)
            csa = aux2(5,iadj,3)
        endif
    
        ! Compute the P-wave parts (a1 downward, a2 upward)
        do j = 1, meqn
            wave(j,1) = 0.d0
            wave(j,2) = 0.d0
        end do
        s(1) = -cpb
        s(2) = cpa

        a1 = (cp*dsig_ii + bulk*du_i) / (bulk*cpb + bulkb*cp)
        a2 = (cp*dsig_ii - bulk*du_i) / (bulk*cpa + bulka*cp)

        wave(sig_ii,1) = a1 * bulkb
        wave(sig_kk1,1) = a1 * lamb
        wave(sig_kk2,1) = a1 * lamb
        wave(u_i,1) = a1 * cpb

        wave(sig_ii,2) = a2 * bulka
        wave(sig_kk1,2) = a2 * lama
        wave(sig_kk2,2) = a2 * lama
        wave(u_i,2) = -a2 * cpa

        ! Compute the S-wave parts (a3,a4 downward, a5,a6 upward)
        do j = 1, meqn
            wave(j,3) = 0.d0
            wave(j,4) = 0.d0
            wave(j,5) = 0.d0
            wave(j,6) = 0.d0
        end do
        s(3) = -csb
        s(4) = -csb
        s(5) = csa
        s(6) = csa

        det = mub*cs + mu*csb
        if (det > 1.e-10) then
            a3 = (cs*dsig_ij1 + mu*du_j1) / det
            a4 = (cs*dsig_ij2 + mu*du_j2) / det

            wave(sig_ij1,3) = a3 * mub
            wave(u_j1,3) = a3 * csb

            wave(sig_ij2,4) = a4 * mub
            wave(u_j2,4) = a4 * csb
        end if

        det = mua*cs + mu*csa
        if (det > 1.e-10) then
            a5 = (cs*dsig_ij1 - mu*du_j1) / det
            a6 = (cs*dsig_ij2 - mu*du_j2) / det

            wave(sig_ij1,5) = a5 * mua
            wave(u_j1,5) = -a5 * csa

            wave(sig_ij2,6) = a6 * mua
            wave(u_j2,6) = -a6 * csa
        end if



!        # Compute downward and upward flux difference:
!       Remember that s1,s3,s4 < 0 and s2,s5,s6 > 0
    
        do j=1,meqn
            bmasdq(j,i) = s(1)*wave(j,1) + s(3)*wave(j,3) + s(4)*wave(j,4)
            bpasdq(j,i) = s(2)*wave(j,2) + s(5)*wave(j,5) + s(6)*wave(j,6)
        end do

    end do

    return
end subroutine rpt3


