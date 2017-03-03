subroutine rp1(maxm,num_eqn,num_waves,num_aux,num_ghost,num_cells,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

! Riemann solver for the 1D shallow water equations:

! (h)_t + (h*u)_x = 0
! (hu)_t + (h*u^2 + 1/2*grav*h^2)_x = -grav*h*(b)_x
!
! using the f-wave algorithm and Roe's approximate Riemann solver.

! waves:     2
! equations: 2

! Conserved quantities:
!       1 depth
!       2 momentum

! Auxiliary fields:
!       1 bathymetry

! The gravitational constant grav should be in the common block cparam.
!
! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

    implicit none 
    
    integer, intent(in) :: maxm, num_eqn, num_waves, num_ghost, num_aux, num_cells
    real(kind=8), intent(in) :: ql(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: qr(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: amdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: apdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    
    real(kind=8) :: grav
    common /cparam/ grav
    
    real(kind=8) :: hl, ul, hr, ur, hbar, uhat, chat, bl, br
    real(kind=8) :: R(2,2), L(2,2)
    real(kind=8) :: fluxdiff(2), beta(2)
    
    integer :: i, j, k
    
    
    do i=2-num_ghost,num_cells+num_ghost
        ! # Left states
        hl = qr(1,i-1)
        ul = qr(2,i-1)/hl
        bl = auxr(1,i-1)
        
        ! # Right states
        hr = ql(1,i)
        ur = ql(2,i)/hr
        br = auxl(1,i)
        
        ! # Roe average states
        hbar = 0.5*(hr+hl)
        uhat = (dsqrt(hl)*ul + dsqrt(hr)*ur)/(dsqrt(hl)+dsqrt(hr))
        chat = dsqrt(grav*hbar)
        
        ! # Flux differences
        fluxdiff(1) = (hr*ur) - (hl*ul)
        fluxdiff(2) = (hr*ur*ur + 0.5*grav*hr**2) - (hl*ul*ul + 0.5*grav*hl**2) + grav*hbar*(br-bl)
        
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
        do j=1,num_waves
            beta(j) = 0.d0
            do k=1,num_eqn
                beta(j) = beta(j) + L(j,k)*fluxdiff(k)
            enddo
        enddo

        ! # Flux waves
        do j=1,num_eqn
            do k=1,num_waves
                fwave(j,k,i) = beta(k)*R(j,k)
            enddo
        enddo
        
        
        ! # Fluctuations
        do j=1,num_eqn
            amdq(j,i) = 0.d0
            apdq(j,i) = 0.d0
            do k=1,num_eqn
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
