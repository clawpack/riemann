! =====================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================

! f-wave Riemann solver for the p-system of Lagrangian compressible gas dynamics in 1d
! Here we assume that the entropy may be spatially varying but constant in time.
! The approach taken is essentially identical to that described in
!
! R. J. LeVeque, 2002, Finite-volume methods for non-linear elasticity in heterogeneous media
!
! The system here is essentially a special case of the one discussed there, with
! \rho(x)=1 and a different equation of state.

!   v_t - u_x = 0
!   u_t + p(v,s(x))_x =0
! where v = specific volume, u = velocity, p=pressure, s=entropy

! waves: 2
! equations: 2
! aux fields: 1

! Conserved quantities:
!       1 specific volume
!       2 velocity

! Auxiliary variables:
!       1 entropy

! function p(v,i) gives pressure in ith cell
! function pp(v,i) gives d/(d v) of pressur

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

    integer, intent(in) :: maxm, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: fwave
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s

    double precision :: enti, entim, Vi, Vim, ui, uim, pi, pim, ci, cim, zi, zim
    double precision :: dp, du, b1, b2
    double precision :: pstar, vstar, gamma
    double precision :: dpdVi, dpdVim
    integer i, m

    common /cparam/  pstar, vstar, gamma

!     # split the jump in q at each interface into waves

    do i = 2-mbc, mx+mbc
        enti = auxl(1,i)
        entim = auxr(1,i-1)
        Vi = ql(1,i)
        Vim = qr(1,i-1)
        ui = ql(2,i)
        uim = qr(2,i-1)

    !        #linearize on each side:

        dpdVi  = -gamma*pstar*vstar*dexp(enti)/(Vi**(gamma+1))
        dpdVim = -gamma*pstar*vstar*dexp(entim)/(Vim**(gamma+1))

        ci = dsqrt(-dpdVi)
        cim = dsqrt(-dpdVim)
        zi = ci
        zim = cim

        du = ui - uim
        pi = pstar*dexp(enti) * (vstar/Vi)**gamma 
        pim = pstar*dexp(entim) * (vstar/Vim)**gamma 
        dp = pi - pim
        b1 = ( dp-zi*du ) / (zim + zi)
        b2 = -(zim*du + dp) / (zim + zi)
    
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
