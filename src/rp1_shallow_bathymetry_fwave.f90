! ================================================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ================================================================================

! # Solve Riemann problems for the 1D shallow water equations
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x = -grav*h*(b)_x
! #
! # using f-wave algorithm and Roe's approximate Riemann solver.  
! #
! # With the f-wave approach the source term in the discharge equation 
! # should be treated here. However, its contribution has been taken into account 
! # in the tfluctf-wave.f90 subroutine, which solves a Riemann problem at the interface.!
! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
! # On output, wave contains the waves, 
! #            s the speeds, 
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
! # Note that the i'th Riemann problem has left state qr(:,i-1)
! #                                    and right state ql(:,i)
! # From the basic clawpack routine step1, rp is called with ql = qr = q.


    implicit none 
    
    integer :: maxm, meqn, mwaves, mbc, maux, mx
    double precision :: ql(meqn, 1-mbc:maxm+mbc)
    double precision :: qr(meqn, 1-mbc:maxm+mbc)
    double precision :: auxl(maux, 1-mbc:maxm+mbc)
    double precision :: auxr(maux, 1-mbc:maxm+mbc)
    double precision :: s(mwaves, 1-mbc:maxm+mbc)
    double precision :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    double precision :: amdq(meqn, 1-mbc:maxm+mbc)
    double precision :: apdq(meqn, 1-mbc:maxm+mbc)
    
    
    double precision :: grav
    common /cparam/ grav
    
    double precision :: hl, ul, hr, ur, hbar, uhat, chat, bl, br
    double precision :: R(2,2)
    double precision :: L(2,2)
    double precision :: fluxDiff(2), beta(2)
    
    integer :: i, j, k
    
    
    do i=2-mbc,mx+mbc
        ! # Left states
        hl = qr(1,i-1)
        ul = qr(2,i-1)/hl
        bl = auxr(1,i-1)
        
        ! # Right states
        hr = ql(1,i)
        ur = ql(2,i)/hr
        br = auxl(1,i)
        
        ! # Average states (they come from the Roe's linearization)
        hbar = 0.5*(hr+hl)
        uhat = (dsqrt(hl)*ul + dsqrt(hr)*ur)/(dsqrt(hl)+dsqrt(hr))
        chat = dsqrt(grav*hbar)
        
        ! # Flux differences
        fluxDiff(1) = (hr*ur) - (hl*ul)
        fluxDiff(2) = (hr*ur*ur + 0.5*grav*hr**2) - (hl*ul*ul + 0.5*grav*hl**2) + grav*hbar*(br-bl)
        
        ! # Wave speeds
        s(1,i) = uhat-chat
        s(2,i) = uhat+chat
        
        ! # Right eigenvectors (column)
        R(1,1) = 1.d0
        R(2,1) = uhat-chat
        
        R(1,2) = 1.d0
        R(2,2) = uhat+chat
    
        ! # Left eigenvectors (rows)
        L(1,1) = (chat+uhat)/(2.d0*chat)
        L(2,1) = (chat-uhat)/(2.d0*chat)
        
        L(1,2) = -1d0/(2.d0*chat)
        L(2,2) = 1.d0/(2.d0*chat)
        
        ! # Coefficients beta
        do j=1,mwaves
            beta(j) = 0.d0
            do k=1,meqn
                beta(j) = beta(j) + L(j,k)*fluxDiff(k)
            enddo
        enddo

        ! # Flux waves
        do j=1,meqn
            do k=1,mwaves
                fwave(j,k,i) = beta(k)*R(j,k)
            enddo
        enddo
        
        
        ! # Fluctuations
        do j=1,meqn
            amdq(j,i) = 0.d0
            apdq(j,i) = 0.d0
            do k=1,meqn
                if (s(k,i) .lt. 0.d0) then
                    amdq(j,i) = amdq(j,i) + fwave(j,k,i)
                elseif (s(k,i) .ge. 0d0) then
                    apdq(j,i) = apdq(j,i) + fwave(j,k,i)
                else
                    amdq(j,i) = amdq(j,i) + 0.5*fwave(j,k,i)
                    apdq(j,i) = apdq(j,i) + 0.5*fwave(j,k,i)
                endif
            enddo
        enddo
    enddo
    
    return
end subroutine rp1
