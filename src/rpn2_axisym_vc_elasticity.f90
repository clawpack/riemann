subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

    ! Riemann solver for the elasticity equations in 2d (axisymmetric), with varying
    ! material properties rho, lambda, and mu 
    !
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
    ! Contents of ql and qr:
    ! 
    ! q(1,:) = sigma_xx
    ! q(2,:) = sigma_yy
    ! q(3,:) = sigma_xy
    ! q(4,:) = sigma_tt
    ! q(5,:) = u (radial)
    ! q(6,:) = v (orthogonal)
    ! 
    ! auxl and auxr hold corresponding slice of the aux array:
    ! Here it is assumed that auxl=auxr gives the cell values
    ! for this slice.
    ! 
    !  auxl(1,i) = rho, density
    !  auxl(2,i) = lambda 
    !  auxl(3,i) = mu
    !  auxl(4,i) = cp, P-wave speed 
    !  auxl(5,i) = cs, S-wave speed 
    !
    !
    ! On output, wave contains the waves,
    !            s the speeds,
    !            amdq the  left-going flux difference  A^- \Delta q
    !            apdq the right-going flux difference  A^+ \Delta q
    !
    ! Note that the waves are *not* in order of increasing lambda.
    ! Instead the 1- and 2-waves are the P-waves and the 3- and 4-waves
    ! are the S-waves.   (The 5th wave has speed zero and is not used.)

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

    integer :: sig_ii, sig_kk, sig_ij, u_i, u_j, m, i, j
    double precision :: dsig_ii, dsig_ij, du_i, du_j
    double precision :: lamr, mur, bulkr, cpr, csr
    double precision :: laml, mul, bulkl, cpl, csl
    double precision :: det, a1, a2, a3, a4

    ! set sig_ii and u_i to correspond to the normal direction
    ! set sig_kk to the diagonal element of stress in the orthogonal direction
    ! set sig_ij to the off-diagonal element of stress
    ! set u_j to the orthoginal direction

    if (ixy.eq.1) then
        sig_ii = 1
        sig_kk = 2
        sig_ij = 3
        u_i = 5
        u_j = 6
    else
        sig_ii = 2
        sig_kk = 1
        sig_ij = 3
        u_i = 6
        u_j = 5
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
        dsig_ii = ql(sig_ii,i) - qr(sig_ii,i-1)
        dsig_ij = ql(sig_ij,i) - qr(sig_ij,i-1)
        du_i = ql(u_i,i) - qr(u_i,i-1)
        du_j = ql(u_j,i) - qr(u_j,i-1)

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

        ! P-wave strengths:
        det = bulkl*cpr + bulkr*cpl
        if (det.eq.0.d0) then
            write(6,*) 'det=0 in rpn2'
            stop 
        endif
        a1 = (cpl*dsig_ii - bulkl*du_i) / det
        a2 = (cpr*dsig_ii + bulkr*du_i) / det

        ! S-wave strengths:
        det = mul*csr + mur*csl
        if (det.eq.0.d0) then
            ! no s-waves
            a3 = 0.d0
            a4 = 0.d0
        else
            a3 = (csl*dsig_ij - mul*du_j) / det
            a4 = (csr*dsig_ij + mur*du_j) / det
        endif

        ! 5th wave has velocity 0 so is not computed or propagated.


        ! Compute the P-waves.
        wave(sig_ii,1,i) = a1 * bulkr
        wave(sig_kk,1,i) = a1 * lamr
        wave(sig_ij,1,i) = 0.d0
        wave(4,1,i) = a1 * lamr
        wave(u_i,1,i) = -a1 * cpr
        wave(u_j,1,i) = 0.d0
        s(1,i) = cpr

        wave(sig_ii,2,i) = a2 * bulkl
        wave(sig_kk,2,i) = a2 * laml
        wave(sig_ij,2,i) = 0.d0
        wave(4,2,i) = a2 * laml
        wave(u_i,2,i) = a2 * cpl
        wave(u_j,2,i) = 0.d0
        s(2,i) = -cpl

        ! Compute the S-waves
        wave(sig_ii,3,i) = 0.d0
        wave(sig_kk,3,i) = 0.d0
        wave(sig_ij,3,i) = a3 * mur
        wave(4,3,i) = 0.d0
        wave(u_i,3,i) = 0.d0
        wave(u_j,3,i) = -a3 * csr
        s(3,i) = csr

        wave(sig_ii,4,i) = 0.d0
        wave(sig_kk,4,i) = 0.d0
        wave(sig_ij,4,i) = a4 * mul
        wave(4,4,i) = 0.d0
        wave(u_i,4,i) = 0.d0
        wave(u_j,4,i) = a4 * csl
        s(4,i) = -csl

        ! compute the leftgoing and rightgoing flux differences:
        ! Note s(i,1),s(i,3) > 0   and   s(i,2),s(i,4) < 0.
        do m=1,meqn
            amdq(m,i) = s(2,i)*wave(m,2,i) + s(4,i)*wave(m,4,i)
            apdq(m,i) = s(1,i)*wave(m,1,i) + s(3,i)*wave(m,3,i)
        enddo
    enddo

    return
end subroutine rpn2
