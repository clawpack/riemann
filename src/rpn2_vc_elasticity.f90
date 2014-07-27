subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

! Riemann solver for the elasticity equations in 2d, with varying
! material properties rho, lambda, and mu 
!
! waves: 4
! equations: 5
! aux fields: 5

! Conserved quantities:
!       1 sigma_11
!       2 sigma_22
!       3 sigma_12
!       4 u
!       5 v

! Auxiliary variables:
!       1  density
!       2  lamda
!       3  mu
!       4  cp
!       5  cs

! Note that although there are 5 eigenvectors, one eigenvalue
! is always zero and so we only need to compute 4 waves.	
! 
! solve Riemann problems along one slice of data.
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr
!
! This data is along a slice in the x-direction if ixy=1 
!                            or the y-direction if ixy=2.
!
! Here it is assumed that auxl=auxr gives the cell values
! for this slice.
!
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!
! Note that the waves are *not* in order of increasing lambda.
! Instead the 1- and 2-waves are the P-waves and the 3- and 4-waves
! are the S-waves.   (The 5th wave would have speed zero and is not computed.)

    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx, maux
    double precision, intent(in) :: ql, qr, auxl, auxr
    double precision, intent(out) :: wave, s, amdq, apdq

    dimension wave( meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn, 1-mbc:maxm+mbc)
    dimension amdq(meqn, 1-mbc:maxm+mbc)
    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)

    integer :: ksig11, ksig22, ku, kv, i, m
    double precision :: dsig11, dsig22, dsig12, du, dv
    double precision :: alamr, amur, bulkr, cpr, csr
    double precision :: alaml, amul, bulkl, cpl, csl
    double precision :: det, a1, a2, a3, a4

    ! set ku to point to  the component of the system that corresponds
    ! to velocity in the direction of this slice, kv to the orthogonal
    ! velocity.  Similarly ksig11 and ksig22 point to normal stresses.
    ! 3rd component is always shear stress sig12.

    if (ixy.eq.1) then
        ksig11 = 1
        ksig22 = 2
        ku = 4
        kv = 5
    else
        ksig11 = 2
        ksig22 = 1
        ku = 5
        kv = 4
    endif

    ! note that notation for u and v reflects assumption that the 
    ! Riemann problems are in the x-direction with u in the normal
    ! direciton and v in the orthogonal direcion, but with the above
    ! definitions of ku and kv the routine also works with ixy=2

    ! split the jump in q at each interface into waves
    ! The jump is split into leftgoing waves traveling at speeds -cp, -cs
    ! relative to the material properties to the left of the interface,
    ! and rightgoing waves traveling at speeds +cp, +cs
    ! relative to the material properties to the right of the interface,

    do i = 2-mbc, mx+mbc
        dsig11 = ql(ksig11,i) - qr(ksig11,i-1)
        dsig22 = ql(ksig22,i) - qr(ksig22,i-1)
        dsig12 = ql(3,i) - qr(3,i-1)
        du = ql(ku,i) - qr(ku,i-1)
        dv = ql(kv,i) - qr(kv,i-1)

        ! material properties in cells i (on right) and i-1 (on left):

        alamr = auxl(2,i)
        amur = auxl(3,i)
        bulkr = alamr + 2.d0*amur
        cpr = auxl(4,i)
        csr = auxl(5,i)

        alaml = auxr(2,i-1)
        amul = auxr(3,i-1)
        bulkl = alaml + 2.d0*amul
        cpl = auxr(4,i-1)
        csl = auxr(5,i-1)

        ! P-wave strengths:
        det = bulkl*cpr + bulkr*cpl
        if (det.eq.0.d0) then
            write(6,*) 'det=0 in rpn2'
            stop 
        endif
        a1 = (cpr*dsig11 + bulkr*du) / det
        a2 = (cpl*dsig11 - bulkl*du) / det

        ! S-wave strengths:
        det = amul*csr + amur*csl
        if (det.eq.0.d0) then
            ! no s-waves
            a3 = 0.d0
            a4 = 0.d0
        else
            a3 = (csr*dsig12 + amur*dv) / det
            a4 = (csl*dsig12 - amul*dv) / det
        endif

        ! 5th wave has velocity 0 so is not computed or propagated.


        ! Compute the waves.

        wave(ksig11,1,i) = a1 * bulkl
        wave(ksig22,1,i) = a1 * alaml
        wave(3,1,i)  = 0.d0
        wave(ku,1,i) = a1 * cpl
        wave(kv,1,i) = 0.d0
        s(1,i) = -cpl

        wave(ksig11,2,i) = a2 * bulkr
        wave(ksig22,2,i) = a2 * alamr
        wave(3,2,i)  = 0.d0
        wave(ku,2,i) = -a2 * cpr
        wave(kv,2,i) = 0.d0
        s(2,i) = cpr

        wave(ksig11,3,i) = 0.d0
        wave(ksig22,3,i) = 0.d0
        wave(3,3,i)  = a3*amul
        wave(ku,3,i) = 0.d0
        wave(kv,3,i) = a3*csl
        s(3,i) = -csl

        wave(ksig11,4,i) = 0.d0
        wave(ksig22,4,i) = 0.d0
        wave(3,4,i)  = a4*amur
        wave(ku,4,i) = 0.d0
        wave(kv,4,i) = -a4*csr
        s(4,i) = csr


        ! compute the leftgoing and rightgoing flux differences:
        ! Note s(i,1),s(i,3) < 0   and   s(i,2),s(i,4) > 0.
        do m=1,meqn
            amdq(m,i) = s(1,i)*wave(m,1,i) + s(3,i)*wave(m,3,i)
            apdq(m,i) = s(2,i)*wave(m,2,i) + s(4,i)*wave(m,4,i)
        enddo
    enddo

    return
end subroutine rpn2
