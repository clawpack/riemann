subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ============================================================================
!  Solves normal Riemann problem for the multilayer shallow water equations in
!  2D with topography:
!    (rho_1 h_1)_t + (rho_1 h_1 u_1)_x + (rho_1 h_1 v_1)_y = 0
!    (rho_1 h_1 u_1)_t + (rho_1 h_1 u_1^2 + 1/2 g rho_1 h_1^2)_x + (rho_1 h_1 u_1 v_1)_y = -g rho_1 h_1(r(h_2)_x + B_x)
!    (rho_1 h_1 v_1)_t + (rho_1 h_1 u_1 v_1)_x + (rho_1 h_1 v_1^2 + 1/2 g rho_1 h_1^2)_y = -g rho_1 h_1(r(h_2)_y + B_y)
!    (rho_2 h_2)_t + (rho_2 h_2 u_2)_x + (rho_2 h_2 v_2)_y = 0
!    (rho_2 h_2 u_2)_t + (rho_2 h_2 u_2^2 + 1/2 g rho_2 h_2^2)_x + (rho_2 h_2 u_2 v_2)_y = -g rho_2 h_2(h_1 + B)_x
!    (rho_2 h_2 v_2)_t + (rho_2 h_2 u_2 v_2)_x + (rho_2 h_2 v_2^2 + 1/2 g rho_2 h_2^2)_y = -g rho_2 h_2(h_1 + B)_y
!
!  On input, ql contains the state vector at the left edge of each cell and qr
!  contains the state vector at the right edge of each cell
!
!           |            |          
!    qr(i-1)|ql(i)  qr(i)|ql(i+1)   
! ----------|------------|----------
!    i-1          i         i+1
!
!  The i-1/2 Riemann problem has left state qr(i-1) and right state ql(i)
!
!  If ixy == 1 then the sweep direction is x, ixy == 2 implies the y direction
! 
!  Wet-Dry Interface Cases handled
!    Dry State Table    #  Handled?    Line #     Comment
!  L(1) L(2) R(1) R(2)
!   F    F    F    F    4  [X]        XXX        Full two-layer case
!
!   T    F    F    F    1  [ ]                   Dry state top layer
!   F    T    F    F    1  [X]        XXX        Lower left dry state
!   F    F    T    F    1  [ ]                   Dry state top layer
!   F    F    F    T    1  [X]        XXX        Lower right dry state
!
!   T    T    F    F    2  [X]        XXX        Left state completely dry
!   F    T    T    F    2  [ ]
!   F    F    T    T    2  [X]        240        Right state completely dry
!   T    F    F    T    2  [ ]
!
!   T    F    T    F    2  [X]        182        Single layer all-wet bottom
!   F    T    F    T    2  [X]        182        Single layer all-wet top
!
!   T    T    T    F    3  [X]        195        Single layer dry-state
!   F    T    T    T    3  [X]        195        Single layer dry-state
!   T    F    T    T    3  [X]        195        Single layer dry-state
!   T    T    F    T    3  [X]        195        Single layer dry-state
!
!   T    T    T    T    4  [X]        173        All dry
! ============================================================================

    use amr_module, only: mcapa

    use geoclaw_module, only: g => grav, rho, pi, earth_radius

    use multilayer_module, only: num_layers, eigen_func
    use multilayer_module, only: dry_tolerance, aux_layer_index, r
    use multilayer_module, only: eigen_method, inundation_method
    use multilayer_module, only: eigen_func, inundation_eigen_func

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx,maux
    real(kind=8), dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    real(kind=8), dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr

    ! Output arguments
    real(kind=8), dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: fwave
    real(kind=8), dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    real(kind=8), dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: apdq, amdq

    ! Counters
    integer :: i, j, mw, info
    integer :: n_index, t_index, layer_index
    
    ! Physics
    real(kind=8) :: dxdc
    
    ! State variables
    real(kind=8), dimension(num_layers) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
    real(kind=8), dimension(num_layers) :: u_l, u_r, v_l, v_r
    real(kind=8), dimension(num_layers) :: h_ave, momentum_transfer
    real(kind=8), dimension(num_layers) :: h_hat_l, h_hat_r
    real(kind=8) :: b_l, b_r, flux_transfer_l, flux_transfer_r, lambda(6)

    ! real(kind=8) :: advected_speed, eta_l, eta_r, gamma_l, gamma_r, kappa_l, kappa_r, w_normal, w_transverse

    ! Solver variables
    integer :: num_dry_states
    real(kind=8), dimension(num_layers) :: eigen_h_l, eigen_h_r
    real(kind=8), dimension(num_layers) :: eigen_u_l, eigen_u_r
    real(kind=8), dimension(num_layers) :: eigen_v_l, eigen_v_r
    real(kind=8), dimension(num_layers) :: flux_h_l, flux_h_r
    real(kind=8), dimension(num_layers) :: flux_hu_l, flux_hu_r
    real(kind=8), dimension(num_layers) :: flux_hv_l, flux_hv_r
    real(kind=8), dimension(num_layers) :: flux_u_l, flux_u_r
    real(kind=8), dimension(num_layers) :: flux_v_l, flux_v_r
    real(kind=8), dimension(6) :: delta, flux_r, flux_l, pivot
    real(kind=8), dimension(6,6) :: eig_vec, A
    real(kind=8) :: beta(6), alpha(4), fw(3, 3), sw(3)
    logical, dimension(num_layers) :: dry_state_l, dry_state_r
    logical :: inundation

    external dgesv

    interface
        subroutine solve_single_layer_rp(layer_index, h_l, h_r, hu_l, hu_r, hv_l, hv_r, b_l, b_r, fw, sw)
            use geoclaw_module, only: g => grav
            use multilayer_module, only: dry_tolerance
            implicit none
            ! Input
            integer, intent(in) :: layer_index
            real(kind=8), intent(in), dimension(2) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
            real(kind=8), intent(in) :: b_l, b_r
            ! Output
            real(kind=8), intent(in out) :: fw(3, 3), sw(3)
        end subroutine solve_single_layer_rp
    end interface
    
    ! Initialize output variables
    amdq = 0.d0
    apdq = 0.d0
    
    ! Set normal direction
    if (ixy == 1) then
        n_index = 2
        t_index = 3
    else
        n_index = 3
        t_index = 2
    endif

    ! ========================================================================
    ! Loop through Riemann problems
    ! ========================================================================
    do i=2-mbc,mx+mbc
        dry_state_l = .false.
        dry_state_r = .false.
        inundation = .false.
        
        ! Parse states and set appropriate zeros
        ! Note that the "u-direction" is the direction of sweeping which 
        ! could actually be the x or y-directions depending on ixy
        
        do j=1,2
            layer_index = 3*(j-1)
            h_l(j) = qr(layer_index+1,i-1) / rho(j)
            hu_l(j) = qr(layer_index+n_index,i-1) / rho(j)
            hv_l(j) = qr(layer_index+t_index,i-1) / rho(j)
            
            h_r(j) = ql(layer_index+1,i) / rho(j)
            hu_r(j) = ql(layer_index+n_index,i) / rho(j)
            hv_r(j) = ql(layer_index+t_index,i) / rho(j)
            
            h_hat_l(j) = auxr(j+aux_layer_index-1,i-1)
            h_hat_r(j) = auxl(j+aux_layer_index-1,i)
            
            h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))

            ! Check for dry states
            if (h_l(j) < dry_tolerance(j)) then
                dry_state_l(j) = .true.
                hu_l(j) = 0.d0
                hv_l(j) = 0.d0
                u_l(j) = 0.d0
                v_l(j) = 0.d0
            else
                u_l(j) = hu_l(j) / h_l(j)
                v_l(j) = hv_l(j) / h_l(j)
            endif
            if (h_r(j) < dry_tolerance(j)) then
                dry_state_r(j) = .true.
                hu_r(j) = 0.d0
                hv_r(j) = 0.d0
                u_r(j) = 0.d0
                v_r(j) = 0.d0
            else
                u_r(j) = hu_r(j) / h_r(j)
                v_r(j) = hv_r(j) / h_r(j)
            endif
        enddo
        
        b_l = auxr(1,i-1)
        b_r = auxl(1,i)

        ! For ease of checking below, count up number of dry states
        num_dry_states = 0
        do mw=1,2
            if (dry_state_l(mw)) then
                num_dry_states = num_dry_states + 1
            end if
            if (dry_state_r(mw)) then
                num_dry_states = num_dry_states + 1
            end if
        end do

        ! ===============================
        !  Completely dry cell - (T T T T)
        ! ===============================
        if (num_dry_states == 4) then
            s(:,i) = 0.d0
            fwave(:,:,i) = 0.d0

        ! ===============================================
        !  Single-layer problem - 
        !   (T F T F) or (F T F T) or (F T T T) or
        !   (T F T T) or (T T F T) or (T T T F)
        ! ===============================================
        !  Here check explicitly for the completely wet single-layer cases and 
        !  rely on the count to determine the other cases
        else if ((      dry_state_l(1) .and. .not. dry_state_l(2) .and.        &  ! T F T F
                        dry_state_r(1) .and. .not. dry_state_r(2))       .or.  &
                 (.not. dry_state_l(1) .and.       dry_state_l(2) .and.        &  ! F T F T
                  .not. dry_state_r(1) .and.       dry_state_r(2))       .or.  &
                 num_dry_states == 3) then

            ! Set wet layer index so that we can handle both the bottom and top
            ! layers being dry while the other is wet
            if (.not.dry_state_l(1) .or. .not.dry_state_r(1)) then
                layer_index = 1
            else if (.not.dry_state_l(2) .or. .not. dry_state_r(2)) then
                layer_index = 2
            else
                print *, "Invalid dry layer state reached."
                print *, "dry states: ", dry_state_l, dry_state_r
                print *, "        left            |             right"
                print *, "====================================================="
                print "(2d16.8)", h_l(1), h_r(1)
                print "(2d16.8)", hu_l(1), hu_r(1)
                print "(2d16.8)", hv_l(1), hv_r(1)
                print "(2d16.8)", h_l(2), h_r(2)
                print "(2d16.8)", hu_l(2), hu_r(2)
                print "(2d16.8)", hv_l(2), hv_r(2)
                print "(2d16.8)", b_l, b_r
                stop
            end if

            call solve_single_layer_rp(layer_index, h_l, h_r, hu_l, hu_r,      &
                                                    hv_l, hv_r, b_l, b_r,      &
                                                    fw, sw)

            ! Update speeds and waves
            ! Note that we represent all the waves in the first three arrays
            ! so it does not directly correspond to the two-layer case's wave
            ! structure
            s(:, i) = 0.d0
            fwave(:, :, i) = 0.d0
            
            if (layer_index == 1) then
                s(1:3, i) = sw(:)
                fwave(1, 1:3, i) = fw(1, :) * rho(layer_index)
                fwave(n_index, 1:3, i) = fw(2, :) * rho(layer_index)
                fwave(t_index, 1:3, i) = fw(3, :) * rho(layer_index)
            else
                s(4:6, i) = sw(:)
                fwave(1, 4:6, i) = fw(1, :) * rho(layer_index)
                fwave(n_index, 4:6, i) = fw(2, :) * rho(layer_index)
                fwave(t_index, 4:6, i) = fw(3, :) * rho(layer_index)
            end if

            ! Go on to next cell, lat-long and fluctuation calculations are 
            ! outside of this loop

        ! ======================================================================
        !  Multi-layer system must be solved
        !   In each case a special eigen system is solved and then the flux
        !     difference is evaluated.
        !   Note that the parameter *eigen_method* controls the method being 
        !     used if the cells are completely wet or there exists a wall 
        !     boundary problem in the bottom layer.  Otherwise the parameter 
        !     *inundation_method* is used.
        ! ======================================================================
        else
            ! By default fill in the eigen and flux evaluation states with their
            ! side values
            if (eigen_method == 1) then
                eigen_h_l = h_hat_l
                eigen_h_r = h_hat_r
            else
                eigen_h_l = h_l
                eigen_h_r = h_r
            end if
            eigen_u_l = u_l
            eigen_u_r = u_r
            eigen_v_l = v_l
            eigen_v_r = v_r

            flux_h_l = h_l
            flux_h_r = h_r
            flux_hu_l = hu_l
            flux_hu_r = hu_r
            flux_hv_l = hv_l
            flux_hv_r = hv_r
            flux_u_l = u_l
            flux_u_r = u_r
            flux_v_l = v_l
            flux_v_r = v_r

            ! Also intialize other flux evaluation stuff
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
            momentum_transfer = 0.d0

            ! ==================================================================
            !  Right state is completely dry - (F F T T)
            if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.          &
                      dry_state_r(1) .and.       dry_state_r(2)) then
    
                ! Inundation occurs
                inundation = sum(h_l) + b_l > b_r
                if (inundation) then
                    print *, "Inundation in this case not yet handled."
                    print *, "  dry_state = ", dry_state_l, dry_state_r
                    stop 

                ! Wall boundary
                else
                    ! Wall state - Mirror left state onto right
                    if (eigen_method /= 1) eigen_h_r = h_l
                    eigen_u_r = -u_l
                    eigen_v_r =  v_r

                    ! Flux evaluation
                    flux_h_r = h_l
                    flux_hu_r = -hu_l
                    flux_u_r = -u_l
                    flux_hv_r = hv_l
                    flux_v_r = v_l

                    flux_transfer_r = 0.d0
                    flux_transfer_l = 0.d0
                    momentum_transfer(1) = 0.d0
                    momentum_transfer(2) = 0.d0
                endif

            ! ==================================================================
            !  Left state is completely dry - (T T F F)
            else if (      dry_state_l(1) .and.       dry_state_l(2) .and.     &
                     .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                ! Inundation
                inundation = sum(h_r) + b_r > b_l
                if (inundation) then
                    print *, "Inundation in this case not yet handled."
                    print *, "  dry_state = ", dry_state_l, dry_state_r
                    stop 

                ! Wall
                else
                    if (eigen_method /= 1) eigen_h_l = h_r
                    eigen_u_l = -u_r
                    eigen_v_l = v_l
                
                    ! Flux evaluation
                    flux_h_l = h_r
                    flux_hu_l = -hu_r
                    flux_u_l = -u_r
                    flux_hv_l = hv_r
                    flux_v_l = v_r

                    flux_transfer_r = 0.d0
                    flux_transfer_l = 0.d0
                    momentum_transfer(1) = 0.d0
                    momentum_transfer(2) = 0.d0

                end if

            ! ==================================================================
            !  Right bottom state is dry - (F F F T)
            else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.     &
                     .not. dry_state_r(1) .and.       dry_state_r(2)) then

                ! Inundation
                inundation = h_l(2) + b_l > b_r
                if (inundation) then
                    if (inundation_method == 1 .or.                            &
                        inundation_method == 4) then
                        eigen_h_r = [h_r(1), 0.d0]
                    else if (inundation_method == 2 .or.                       &
                             inundation_method == 3 .or.                       &
                             inundation_method == 5) then
                        eigen_h_r = [h_r(1), 0.d0]
                    end if
    
                    ! Flux evaluation
                    momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
                    momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
                    flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
                    flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

                ! Wall
                else
                    if (eigen_method /= 1) eigen_h_r = [h_r(1), 0.d0]
                    eigen_u_r = [u_r(1), -u_l(2)]
                    eigen_v_r = [v_r(1),  v_l(2)]

                    ! Flux evaluation
                    flux_h_r(2) = h_l(2)
                    flux_hu_r(2) = -hu_l(2)
                    flux_u_r(2) = -u_l(2)
                    flux_hv_r(2) = hv_l(2)
                    flux_v_r(2) = v_l(2)
                
                    flux_transfer_r = 0.d0
                    flux_transfer_l = 0.d0
                    momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - flux_h_l(2) - b_l)
                    momentum_transfer(2) = 0.d0

                end if

            ! ==================================================================
            !  Left bottom state is dry - (F T F F)
            else if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.     &
                     .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                ! Inundation
                inundation = (h_r(2) + b_r > b_l)
                if (inundation) then
                    if (inundation_method == 1 .or. inundation_method == 5) then
                        eigen_h_l = [h_l(1), 0.d0]
                    else if (inundation_method == 2 .or.                       &
                             inundation_method == 3 .or.                       &
                             inundation_method == 4) then
                        eigen_h_l = [h_l(1), dry_tolerance(1)]
                    end if
    
                    ! Flux evaluation
                    momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
                    momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
                    flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
                    flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

                ! Wall
                else
                    if (eigen_method /= 1) eigen_h_l = [h_l(1), 0.d0]
                    eigen_u_l = [u_l(1), -u_r(2)]
                    eigen_v_l = [v_l(1),  v_r(2)]
    
                    ! Flux evaluation
                    flux_h_l(2) = h_r(2)
                    flux_hu_l(2) = -hu_r(2)
                    flux_u_l(2) = -u_r(2)
                    flux_hv_l(2) = hv_r(2)
                    flux_v_l(2) = v_r(2)
                
                    flux_transfer_r = 0.d0
                    flux_transfer_l = 0.d0
                    momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + flux_h_r(2) - b_l)
                    momentum_transfer(2) = 0.d0

                end if

            ! ==================================================================
            !  All-states are wet - (F F F F)
!             else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.     &
!                      .not. dry_state_r(1) .and. .not. dry_state_r(2)) then
            else if (num_dry_states == 0) then

                ! Nothing to do for eigenspace evaluation

                ! Flux evaulation
                momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
                momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
                flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
                flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

            ! ==================================================================
            !  We do not yet handle this case - F F F F and F F F F 
            else
                print *, "Unhandled dry-state condition reached."
                print *, "dry states: ", dry_state_l, dry_state_r
                print *, "        left            |             right"
                print *, "====================================================="
                print "(2d16.8)", h_l(1), h_r(1)
                print "(2d16.8)", hu_l(1), hu_r(1)
                print "(2d16.8)", hv_l(1), hv_r(1)
                print "(2d16.8)", h_l(2), h_r(2)
                print "(2d16.8)", hu_l(2), hu_r(2)
                print "(2d16.8)", hv_l(2), hv_r(2)
                print "(2d16.8)", b_l, b_r
                stop
            end if

            ! ==================================================================
            !  Compute eigen space
            ! ==================================================================
            if (inundation) then
                call inundation_eigen_func(eigen_h_l, eigen_h_r,               &
                                           eigen_u_l, eigen_u_r,               &
                                           eigen_v_l, eigen_v_r,               &
                                           n_index, t_index,                   &
                                           lambda, eig_vec)
                                
                ! Internal wave corrections
                if (inundation_method == 1) then
                    ! Left bottom state dry
                    if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.  &
                        .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                        s(2,i) = u_r(2) - 2.d0 * sqrt(g*(1.d0-r)*h_r(2))
                        alpha(2) = r * g * h_r(2) / ((s(2,i) - u_r(2))**2 - g*h_r(2))
                        eig_vec(1,2) = 1.d0
                        eig_vec(n_index,2) = s(2,i) 
                        eig_vec(t_index,2) = v_l(1)
                        eig_vec(4,2) = alpha(2)
                        eig_vec(n_index,2) = alpha(2)*s(2,i)
                        eig_vec(t_index,2) = alpha(2)*v_l(2)
                    ! Right bottom state dry
                    else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and. &
                             .not. dry_state_r(1) .and.       dry_state_r(2)) then

                        s(5,i) = u_l(2) + 2.d0 * sqrt(g*(1.d0-r)*h_l(2))
                        alpha(3) = r * g * h_l(2) / ((s(5,i) - u_l(2))**2 - g * h_l(2))
                        
                        eig_vec(1,5) = 1.d0
                        eig_vec(n_index,5) = s(5,i)
                        eig_vec(t_index,5) = v_r(1)
                        eig_vec(4,5) = alpha(3)
                        eig_vec(n_index,5) = s(5,i) * alpha(3)
                        eig_vec(t_index,6) = v_r(2) * alpha(3)
                    end if
                end if
                if (inundation_method /= 5) then
                    ! Left bottom state dry
                    if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.  &
                        .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                        s(1,i) = u_l(1) - sqrt(g*h_l(1))
                        eig_vec(1,1) = 1.d0
                        eig_vec(n_index,1) = s(1,i)
                        eig_vec(t_index,1) = v_l(1)
                        eig_vec(4:6,1) = 0.d0

                    ! Right bottom state dry
                    else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and. &
                             .not. dry_state_r(1) .and.       dry_state_r(2)) then

                        s(6,i) = u_r(1) + sqrt(g*h_r(1))
                        eig_vec(1,6) = 1.d0
                        eig_vec(n_index,6) = s(6,i)
                        eig_vec(t_index,6) = v_r(1)
                        eig_vec(4:6,6) = 0.d0

                    end if

                end if
            else
                call eigen_func(eigen_h_l, eigen_h_r,                    &
                                eigen_u_l, eigen_u_r,                    &
                                eigen_v_l, eigen_v_r,                    &
                                n_index, t_index,                        &
                                lambda, eig_vec)

            end if

            s(:,i) = lambda

            ! ======================================================================
            !  Compute flux differences
            ! ======================================================================
            do j=1,2
                layer_index = 3*(j-1)
                flux_r(layer_index+1) = rho(j) * flux_hu_r(j)
                flux_r(layer_index+n_index) = rho(j) * (flux_h_r(j) * flux_u_r(j)**2 + 0.5d0 * g * flux_h_r(j)**2)
                flux_r(layer_index+t_index) = rho(j) * flux_h_r(j) * flux_u_r(j) * flux_v_r(j)
                
                flux_l(layer_index+1) = rho(j) * flux_hu_l(j)
                flux_l(layer_index+n_index) = rho(j) * (flux_h_l(j) * flux_u_l(j)**2 + 0.5d0 * g * flux_h_l(j)**2)
                flux_l(layer_index+t_index) = rho(j) * flux_h_l(j) * flux_u_l(j) * flux_v_l(j)
            enddo
            ! Add extra flux terms
            flux_r(3 + n_index) = flux_r(3 + n_index) + flux_transfer_r
            flux_l(3 + n_index) = flux_l(3 + n_index) + flux_transfer_l
            
            delta = flux_r - flux_l
                
            ! Momentum transfer and bathy terms
            delta(n_index) = delta(n_index) + momentum_transfer(1)
            delta(n_index+3) = delta(n_index+3) + momentum_transfer(2)

            ! ======================================================================
            ! Project jump in fluxes - Use LAPACK's dgesv routine
            !    N - (int) - Number of linear equations (6)
            !    NRHS - (int) - Number of right hand sides (1)
            !    A - (dp(6,6)) - Coefficient matrix, in this case eig_vec
            !    LDA - (int) - Leading dimension of A (6)
            !    IPIV - (int(N)) - Pivot indices
            !    B - (dp(LDB,NRHS)) - RHS of equations (delta)
            !    LDB - (int) - Leading dimension of B (6)
            !    INFO - (int) - Status of result
            !  Note that the solution (betas) are in delta after the call
            ! ======================================================================
            A = eig_vec ! We need to do this as the return matrix is modified and
                        ! we have to use eig_vec again to compute fwaves
            info = 0       
            call dgesv(6,1,A,6,pivot,delta,6,info)
            if (.not.(info == 0)) then
                print *, "dry states: ", dry_state_l, dry_state_r
                print *, "        left            |             right"
                print *, "====================================================="
                print "(2d16.8)", h_l(1), h_r(1)
                print "(2d16.8)", hu_l(1), hu_r(1)
                print "(2d16.8)", hv_l(1), hv_r(1)
                print "(2d16.8)", h_l(2), h_r(2)
                print "(2d16.8)", hu_l(2), hu_r(2)
                print "(2d16.8)", hv_l(2), hv_r(2)
                print "(2d16.8)", b_l, b_r
                print *,""
                print "(a,i2)","In normal solver: ixy=",ixy
                print "(a,i3)","  Error solving R beta = delta,",info
                print "(a,6d16.8)","  Eigenspeeds: ",s(:,i)
                print "(a)","  Eigenvectors:"
                do j=1,6
                    print "(a,6d16.8)","  ",(eig_vec(j,mw),mw=1,6)
                enddo
                stop
            endif
            beta = delta

            ! ======================================================================
            ! Compute fwaves
            forall(mw=1:mwaves)
                fwave(:,mw,i) = eig_vec(:,mw) * beta(mw)
            end forall
        end if

    end do
    ! == End of Riemann Solver Loop per grid cell ==============================

    ! ==========================================================================
    ! Capacity for mapping from latitude longitude to physical space
    if (mcapa > 0) then
        do i=2-mbc,mx+mbc
            if (ixy == 1) then
                dxdc=(earth_radius*pi/180.d0)
            else
                dxdc=earth_radius*cos(auxl(3,i))*pi/180.d0
            endif

            do mw=1,mwaves
                s(mw,i)=dxdc*s(mw,i)
                fwave(:,mw,i)=dxdc*fwave(:,mw,i)
            enddo
        enddo
    endif

    ! ==========================================================================
    !  Compute fluctuations 
    do i=2-mbc,mx+mbc
        do mw=1,mwaves
            if (s(mw,i) > 0.d0) then
                apdq(:,i) = apdq(:,i) + fwave(:,mw,i)
            else
                amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
            endif
            !h_r(1) = ql(1,i) / rho(1)
            !h_l(1) = qr(1,i-1) / rho(1)
            !h_r(2) = ql(4,i) / rho(2)
            !h_l(2) = qr(4,i-1) / rho(2)
            !dry_state_r(2) = h_r(2) < dry_tolerance(2)
            !dry_state_l(2) = h_l(2) < dry_tolerance(2)
            !rare(1) = h_l(2) + b_l > b_r
            !rare(2) = h_r(2) + b_r > b_l
            !if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(.not.rare(1)) &
            ! .or.(dry_state_l(2).and.(.not.dry_state_r(2)).and.(.not.rare(2)))) then
            !    do m=4,6
            !        if (abs(apdq(m,i)) <= 1.d-8) then
            !            print *,"========================"
            !            print *,"Wave ",mw," equation ",m
            !            print *,"s = ",s(mw,i)
            !            print *,"f = ",fwave(m,mw,i)
            !            print *,"amdq = ",(amdq(m,i))
            !            print *,"apdq = ",(apdq(m,i))
            !            stop "Fluctuation apdq non-zero going into a wall, aborting calculation."
            !        endif
            !    enddo
            !    apdq(4:6,i) = 0.d0
            !endif
        enddo
    enddo

end subroutine rpn2


! ==================================
!  Solve the single-layer equations
! ==================================
! subroutine solve_sinlge_layer_rp(layer_index, h_l, h_r,                        &
!                                               hu_l, hu_r,                      &
!                                               hv_l, hv_r,                      &
!                                               u_l, u_r,                        &
!                                               v_l, v_r,                        &
!                                               b_l, b_r,                        &
!                                               fw, sw)

!     use geoclaw_module, only: g => grav
!     use multilayer_module, only: dry_tolerance

!     implicit none

!     ! Input
!     integer, intent(in) :: layer_index
!     real(kind=8), intent(in out), dimension(2) :: h_l, h_r
!     real(kind=8), intent(in out), dimension(2) :: hu_l, hu_r
!     real(kind=8), intent(in out), dimension(2) :: hv_l, hv_r
!     real(kind=8), intent(in out), dimension(2) :: u_r, u_l
!     real(kind=8), intent(in out), dimension(2) :: v_r, v_l
!     real(kind=8), intent(in out) :: b_l, b_r

!     ! Output
!     real(kind=8), intent(out) :: fw(3, 3), sw(3)

!     ! Local storage
!     integer :: m, mw
!     logical :: rare(2)
!     real(kind=8) :: wall(3), h_star, h_star_test, sm(2)
!     real(kind=8) :: phi_l, phi_r, s_l, s_r, u_hat, c_hat, s_roe(2), s_E(2)

!     ! Algorithm parameters
!     integer, parameter :: MAX_ITERATIONS = 1

!     wall = 1.d0
    
!     ! Calculate momentum fluxes
!     phi_l = 0.5d0 * g * h_l(layer_index)**2      &
!                         + h_l(layer_index) * u_l(layer_index)**2
!     phi_r = 0.5d0 * g * h_r(layer_index)**2      &
!                         + h_r(layer_index) * u_r(layer_index)**2
     
!     ! Check for dry state to right
!     if (h_r(layer_index) < dry_tolerance(layer_index)) then
!         call riemanntype(h_l(layer_index), h_l(layer_index),   &
!                          u_l(layer_index),-u_l(layer_index),   &
!                          h_star, sm(1), sm(2), rare(1), rare(2), 1,    &
!                          dry_tolerance(layer_index), g)
!         h_star_test = max(h_l(layer_index), h_star)
!         ! Right state should become ghost values that mirror left for wall problem
!         if (h_star_test + b_l < b_r) then 
!             wall(2:3) = 0.d0
!             h_r(layer_index) = h_l(layer_index)
!             hu_r(layer_index) = -hu_l(layer_index)
!             b_r = b_l
!             phi_r = phi_l
!             u_r(layer_index) = -u_l(layer_index)
!             v_r(layer_index) = v_l(layer_index)
!         else if (h_l(layer_index) + b_l < b_r) then
!             b_r = h_l(layer_index) + b_l
!         endif
!     ! Check for drystate to left, i.e right surface is lower than left topo
!     else if (h_l(layer_index) < dry_tolerance(layer_index)) then 
!         call riemanntype(h_r(layer_index), h_r(layer_index),   &
!                         -u_r(layer_index), u_r(layer_index),   &
!                          h_star, sm(1), sm(2), rare(1), rare(2), 1,    &
!                          dry_tolerance(layer_index), g)
!         h_star_test = max(h_r(layer_index), h_star)
!         ! Left state should become ghost values that mirror right
!         if (h_star_test + b_r < b_l) then  
!            wall(1:2) = 0.d0
!            h_l(layer_index) = h_r(layer_index)
!            hu_l(layer_index) = -hu_r(layer_index)
!            b_l = b_r
!            phi_l = phi_r
!            u_l(layer_index) = -u_r(layer_index)
!            v_l(layer_index) = v_r(layer_index)
!         else if (h_r(layer_index) + b_r < b_l) then
!            b_l = h_r(layer_index) + b_r
!         endif
!     endif

!     ! Determine wave speeds
!     ! 1 wave speed of left state
!     s_l = u_l(layer_index) - sqrt(g * h_l(layer_index))
!     ! 2 wave speed of right state
!     s_r = u_r(layer_index) + sqrt(g * h_r(layer_index))
    
!     ! Roe average
!     u_hat = (sqrt(g * h_l(layer_index)) * u_l(layer_index)     &
!            + sqrt(g * h_r(layer_index)) * u_r(layer_index))    &
!            / (sqrt(g * h_r(layer_index)) + sqrt(g * h_l(layer_index))) 
!     c_hat = sqrt(g * 0.5d0 * (h_r(layer_index) + h_l(layer_index))) 
!     s_roe(1) = u_hat - c_hat ! Roe wave speed 1 wave
!     s_roe(2) = u_hat + c_hat ! Roe wave speed 2 wave
!     s_E(1) = min(s_l, s_roe(1)) ! Eindfeldt speed 1 wave
!     s_E(2) = max(s_r, s_roe(2)) ! Eindfeldt speed 2 wave
    
!     ! Solve Riemann problem
!     call riemann_aug_JCP(MAX_ITERATIONS, 3, 3,                         &
!                          h_l(layer_index), h_r(layer_index),   &
!                          hu_l(layer_index), hu_r(layer_index), &
!                          hv_l(layer_index), hv_r(layer_index), &
!                          b_l, b_r,                                     &
!                          u_l(layer_index), u_r(layer_index),   &
!                          v_l(layer_index), v_r(layer_index),   &
!                          phi_l, phi_r, s_E(1), s_E(2),           &
!                          dry_tolerance(layer_index), g, sw, fw)
    
!     ! Eliminate ghost fluxes for wall
!     do mw=1,3
!         sw(mw) = sw(mw) * wall(mw)
!         do m=1,3
!            fw(m, mw) = fw(m, mw) * wall(mw)
!         enddo
!     enddo

! end subroutine solve_sinlge_layer_rp

subroutine solve_single_layer_rp(layer_index, h_l, h_r, hu_l, hu_r, hv_l, hv_r, b_l, b_r, fw, sw)

    use geoclaw_module, only: g => grav
    use multilayer_module, only: dry_tolerance

    implicit none

    ! Input
    integer, intent(in) :: layer_index
    real(kind=8), intent(in), dimension(2) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
    real(kind=8), intent(in) :: b_l, b_r

    ! Output
    real(kind=8), intent(in out) :: fw(3, 3), sw(3)

    ! Locals
    integer :: mw
    real(kind=8) :: hL, hR, huL, huR, hvL, hvR, uL, uR, vL, vR, bL, bR, pL, pR
    real(kind=8) :: phiL, phiR, wall(3), drytol
    real(kind=8) :: hstar, hstartest, s1m, s2m, rare1, rare2, sL, sR, uhat, chat, sRoe1, sRoe2, sE1, sE2

    ! Parameters (should be anyway)
    integer :: maxiter

    ! Density used in pressure gradient calculation
    ! This needs to be realistic here as compared to air density so we cannot
    ! use what is stored in the multilayer_module unless real densities are 
    ! being used
    ! TODO - fix this limitation 
    real(kind=8), parameter :: rho = 1025.d0

    drytol = dry_tolerance(layer_index)

    hL = h_l(layer_index)
    hR = h_r(layer_index)
    huL = hu_l(layer_index)
    huR = hu_r(layer_index)
    hvL = hv_l(layer_index)
    hvR = hv_r(layer_index)
    bL = b_l
    bR = b_r
    pL = 0.d0
    pR = 0.d0

    ! ========================================
    !  Begin Snipped Code From rpn2_geoclaw.f
    ! ========================================
    !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         endif

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
!                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
!               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

         !determine wave speeds
         sL=uL-sqrt(g*hL) ! 1 wave speed of left state
         sR=uR+sqrt(g*hR) ! 2 wave speed of right state

         uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR, &
                                          vL,vR,phiL,phiR,pL,pR,sE1,sE2,     &
                                          drytol,g,rho,sw,fw)

!         call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
!     &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

!          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
!     &      bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

!        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)

               fw(1,mw)=fw(1,mw)*wall(mw) 
               fw(2,mw)=fw(2,mw)*wall(mw)
               fw(3,mw)=fw(3,mw)*wall(mw)
         enddo

end subroutine solve_single_layer_rp
