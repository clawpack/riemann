! ==================================================================
subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! ==================================================================

! Riemann solver for the elasticity equations in 3d, with varying
! material properties.
!
! waves: 6
! equations: 9
! aux fields: 5
!
! Conserved quantities:
!       1 sigma_xx
!       2 sigma_yy
!       3 sigma_zz
!       4 sigma_xy
!       5 sigma_xz
!       6 sigma_yz
!       7 u
!       8 v
!       9 w
!
! Auxiliary variables:
!       1 rho
!       2 lambda
!       3 mu
!       4 cp
!       5 cs

! Note that although there are 9 eigenvectors, 3 eigenvalues are
! always zero and so we only need to compute 6 waves.

! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

! Note waves 1-2 are the P-waves, waves 3-6 are the S-waves

    implicit none
    integer, intent(in) :: ixyz, maxm,meqn,mwaves,mbc,mx, maux
    double precision, intent(in) ::   ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) ::   qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: auxl(maux,1-mbc:maxm+mbc)
    double precision, intent(in) :: auxr(maux,1-mbc:maxm+mbc)
    double precision, intent(out) :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision, intent(out) ::    s(mwaves,1-mbc:maxm+mbc)
    double precision, intent(out) :: amdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: apdq(meqn,1-mbc:maxm+mbc)
    integer :: i, j, sig_ii, sig_ij1, sig_ij2, sig_kk1, sig_kk2, u_i, u_j1, u_j2
    double precision :: dsig_ii, dsig_ij1, dsig_ij2, du_i, du_j1, du_j2
    double precision :: laml, mul, bulkl, cpl, csl, lamr, mur, bulkr, cpr, csr
    double precision :: det, a1, a2, a3, a4, a5, a6



!     # set sig_ii to the diagonal component of stress in this direction
!     # set sig_ij1 and sig_ij2 to the off-diagonal components of stress in this direction
!     # set u_i to the velocity component in this direction
!     # set u_j1 and u_j2 to the velocity components orthogonal to this direction

    if (ixyz == 1) then
        sig_ii = 1  ! sig_xx
        sig_kk1 = 2 ! sig_yy
        sig_kk2 = 3 ! sig_zz
        sig_ij1 = 4 ! sig_xy
        sig_ij2 = 5 ! sig_xz
        u_i = 7     ! u
        u_j1 = 8    ! v
        u_j2 = 9    ! w
    else if (ixyz == 2) then
        sig_ii = 2  ! sig_yy
        sig_kk1 = 1 ! sig_xx
        sig_kk2 = 3 ! sig_zz
        sig_ij1 = 4 ! sig_xy
        sig_ij2 = 6 ! sig_yz
        u_i = 8     ! v
        u_j1 = 7    ! u
        u_j2 = 9    ! w
    else if (ixyz == 3) then
        sig_ii = 3  ! sig_zz
        sig_kk1 = 1 ! sig_xx
        sig_kk2 = 2 ! sig_yy
        sig_ij1 = 5 ! sig_xz
        sig_ij2 = 6 ! sig_yz
        u_i = 9     ! w
        u_j1 = 7    ! u
        u_j2 = 8    ! v
    endif

!     # split the jump in q at each interface into waves
!     # The jump is split into 1 leftgoing wave traveling at speed -cp
!     # and 2 leftgoing waves traveling at speed -cs
!     # relative to the material properties to the left of the interface,
!     # and 1 rightgoing wave traveling at speed cp
!     # and 2 rightgoing waves traveling at speed cs
!     # relative to the material properties to the right of the interface,



!     # find a1-a6, the coefficients of the 6 eigenvectors:
    do i = 2-mbc, mx+mbc
        dsig_ii = ql(sig_ii,i) - qr(sig_ii,i-1)
        dsig_ij1 = ql(sig_ij1,i) - qr(sig_ij1,i-1)
        dsig_ij2 = ql(sig_ij2,i) - qr(sig_ij2,i-1)
        du_i = ql(u_i,i) - qr(u_i,i-1)
        du_j1 = ql(u_j1,i) - qr(u_j1,i-1)
        du_j2 = ql(u_j2,i) - qr(u_j2,i-1)

        ! material properties in cells i (on right) and i-1 (on left):
        lamr = auxl(2,i)
        mur = auxl(3,i)
        bulkr = lamr + 2.d0*mur
        cpr = auxl(4,i)
        csr = auxl(5,i)

        laml = auxr(2,i-1)
        mul = auxr(3,i-1)
        bulkl = laml + 2.d0*mul
        cpl = auxr(4,i-1)
        csl = auxr(5,i-1)

        ! Compute the P-waves
        do j = 1, meqn
            wave(j,1,i) = 0.d0
            wave(j,2,i) = 0.d0
        end do
        s(1,i) = -cpl
        s(2,i) = cpr

        det = bulkl*cpr + bulkr*cpl
        if (det < 1.e-10) then
            write(6,*) 'det=0 in rpn3'
            write (6,*) 'cpr', cpr, 'bulkl', bulkl, 'cpl', cpl, 'bulkr', bulkr, 'i', i
            write (6,*) 'test', auxl(1,i)

            stop
        else
            a1 = (cpr*dsig_ii + bulkr*du_i) / det
            a2 = (cpl*dsig_ii - bulkl*du_i) / det

            wave(sig_ii,1,i) = a1 * bulkl
            wave(sig_kk1,1,i) = a1 * laml
            wave(sig_kk2,1,i) = a1 * laml
            wave(u_i,1,i) = a1 * cpl

            wave(sig_ii,2,i) = a2 * bulkr
            wave(sig_kk1,2,i) = a2 * lamr
            wave(sig_kk2,2,i) = a2 * lamr
            wave(u_i,2,i) = -a2 * cpr
        end if


        ! Compute the S-waves
        do j = 1, meqn
            wave(j,3,i) = 0.d0
            wave(j,4,i) = 0.d0
            wave(j,5,i) = 0.d0
            wave(j,6,i) = 0.d0
        end do
        s(3,i) = -csl
        s(4,i) = -csl
        s(5,i) = csr
        s(6,i) = csr

        det = mul*csr + mur*csl
        if (det > 1.e-10) then
            a3 = (csr*dsig_ij1 + mur*du_j1) / det
            a4 = (csr*dsig_ij2 + mur*du_j2) / det
            a5 = (csl*dsig_ij1 - mul*du_j1) / det
            a6 = (csl*dsig_ij2 - mul*du_j2) / det

            wave(sig_ij1,3,i) = a3 * mul
            wave(u_j1,3,i) = a3 * csl

            wave(sig_ij2,4,i) = a4 * mul
            wave(u_j2,4,i) = a4 * csl

            wave(sig_ij1,5,i) = a5 * mur
            wave(u_j1,5,i) = -a5 * csr

            wave(sig_ij2,6,i) = a6 * mur
            wave(u_j2,6,i) = -a6 * csr
        end if


        ! compute the leftgoing and rightgoing flux differences:
        ! Note s1,s3,s4 < 0   and   s2,s5,s6 > 0.

        do j=1,meqn
            amdq(j,i) = s(1,i)*wave(j,1,i) + s(3,i)*wave(j,3,i) + s(4,i)*wave(j,4,i)
            apdq(j,i) = s(2,i)*wave(j,2,i) + s(5,i)*wave(j,5,i) + s(6,i)*wave(j,6,i)
        end do
    end do

    return
end subroutine rpn3

