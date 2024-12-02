! =====================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================

! This version uses interleaved arrays.
! This version outputs f-waves.

! Riemann solver for the nonlinear elastic equations in 1d,
!  variable coefficients
!   eps_t - (m/rho(x))_x = 0
!   m_t - sigma(eps,x)_x =0
! where eps=strain, m=rho*u=momentum

! waves: 2
! equations: 2
! aux fields: 2

! Conserved quantities:
!       1 strain
!       2 momentum

! Auxiliary variables:
!       1 density
!       2 bulk_modulus

! function sigma(eps,i) gives stress-strain relation in ith cell
! function sigmap(eps,i) gives d/(d eps) of sigma
!    For linear:   sigma(eps,i) = K_i * eps, and sigmap(eps,i) = K_i

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! On output, fwave contains the waves as jumps in f,
!            s the speeds,
!            amdq = A^- Delta q,
!            apdq = A^+ Delta q,
!                   the decomposition of the flux difference
!                       f(qr(i-1)) - f(ql(i))
!                   into leftgoing and rightgoing parts respectively.
! 

! Note that the ith Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(in) :: ql(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: qr(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxl(maux,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxr(maux,1-mbc:maxmx+mbc)

    real(kind=8), intent(out) :: s(mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: fwave(meqn, mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: amdq(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: apdq(meqn,1-mbc:maxmx+mbc)

    integer :: i, m
    real(kind=8) :: rhoi, rhoim, epsi, epsim, urhoi, urhoim
    real(kind=8) :: bulki, bulkim, ci, cim, zi, zim, du, dsig, b1, b2
    real(kind=8) :: sigma, sigmap
!     # split the jump in q at each interface into waves

    do i = 2-mbc, mx+mbc
        rhoi = auxl(1,i)
        rhoim = auxr(1,i-1)
        epsi = ql(1,i)
        epsim = qr(1,i-1)
        urhoi = ql(2,i)
        urhoim = qr(2,i-1)

    !        #linearize on each side:

        bulki = sigmap(epsi,i,auxl(2,i),auxl(3,i))
        bulkim = sigmap(epsim,i-1,auxr(2,i-1),auxr(3,i-1))
        ci = dsqrt(bulki / rhoi)
        cim = dsqrt(bulkim / rhoim)
        zi = ci*rhoi
        zim = cim*rhoim

        du = urhoi/rhoi - urhoim/rhoim
        dsig = sigma(epsi,i,auxl(2,i),auxl(3,i)) &
        - sigma(epsim,i-1,auxr(2,i-1),auxr(3,i-1))
        b1 = -(zi*du + dsig) / (zim + zi)
        b2 = -(zim*du - dsig) / (zim + zi)
    !         a1 = b1 / (-cim)
    !         a2 = b2 / ci

    !         estarim = epsim + a1
    !         estari = epsi - a2

    !         ui = urhoi / rhoi
    !         uim = urhoim / rhoim

    !         ustar = urhoim/rhoim + (estarim - epsim) * cim
    !               = urhoi/rhoi - (estari - epsi) * ci
    !               = uim - b1
                
    
    !        # Compute the waves.
    
        fwave(1,1,i) = b1
        fwave(2,1,i) = b1 * zim
        s(1,i) = -cim
    
        fwave(1,2,i) = b2
        fwave(2,2,i) = b2*(-zi)
        s(2,i) = ci

    end do


!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

    do m=1,meqn
        do i = 2-mbc, mx+mbc
            amdq(m,i) = fwave(m,1,i)
            apdq(m,i) = fwave(m,2,i)
        end do
    end do

    return
    end subroutine rp1


!     --------------------------------------------
    double precision function sigma(eps,i,coefA,coefB)
!     --------------------------------------------
    implicit none

    real(kind=8) :: eps, coefA, coefB
    integer :: i

!     # stress-strain relation in ith cell


!     # nonlinear in both layers:
!     sigma = (coefB + coefA)*(dexp(eps) - 1.d0)

!     # nonlinear layer 0:
!     sigma = coefB*eps + coefA*(dexp(eps) - 1.d0)

!     # nonlinear layer 1:
!     sigma = coefB*(dexp(eps) - 1.d0) + coefA*eps

    if (coefA > 0.d0) then
    !         # layer A:
    !         sigma = coefA*eps
        sigma = dexp(coefA*eps) - 1.d0
    !          sigma = coefA*eps + 0.5d0*eps**2
    else
    !         # layer B:
        sigma = coefB*eps
    !         sigma = dexp(coefB*eps) - 1.d0
    endif

    return
    END function


!     --------------------------------------------
    double precision function sigmap(eps,i,coefA,coefB)
!     --------------------------------------------
    implicit none

    real(kind=8) :: eps, coefA, coefB
    integer :: i

!     # derivative of stress-strain relation in ith cell

!     # nonlinear in both layers:
!     sigmap = (coefB + coefA)*dexp(eps)

!     # nonlinear layer 0:
!     sigmap = coefB + coefA*dexp(eps)

!     # nonlinear layer 1:
!     sigmap = coefB*dexp(eps) + coefA

    if (coefA > 0.d0) then
    !         # layer A:
    !         sigmap = coefA
        sigmap = coefA*dexp(coefA*eps)
    !         sigmap = coefA + eps
    else
    !         # layer B:
        sigmap = coefB
    !         sigmap = coefB*dexp(coefB*eps)
    endif


    return
    END function