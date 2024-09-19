subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! ============================================================================
!  Solves transverse Riemann problem for the multilayer shallow water 
!  equations in 2D with topography and wind forcing:
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
!  if imp == 1 then the split is in the negative direction, imp == 2 positive
!
!  Kyle T. Mandli (10-13-2010)
! ==============================================================================
!  Note that this version does not work with PyClaw due to missing parameters
!   from the modules.
! ==============================================================================

    use amr_module, only: mcapa
    
    use geoclaw_module, only: g => grav, rho, earth_radius, pi
    use geoclaw_module, only: coordinate_system

    use multilayer_module, only: num_layers, eigen_method, inundation_method
    use multilayer_module, only: eigen_func, inundation_eigen_func
    use multilayer_module, only: dry_tolerance, aux_layer_index

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx,imp,maux
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(inout) :: asdq
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    
    ! Ouput
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq
    
    ! Local storage
    integer :: i,j,m,mw,n_index,t_index,info
    double precision :: dxdcm,dxdcp
    double precision, dimension(3) :: b
    double precision, dimension(6) :: s,delta,pivot
    double precision, dimension(6,6) :: eig_vec,A
    
    ! Single layer solver storage
    double precision :: ql_sl(3),qr_sl(3)
    double precision :: aux1_sl(3,2),aux2_sl(3,2),aux3_sl(3,2)
    double precision :: asdq_sl(3),bmasdq_sl(3),bpasdq_sl(3)
    
    double precision :: ts(6),teig_vec(6,6)
    
    integer :: layer_index
    logical :: dry_state_l(2),dry_state_r(2)
    double precision :: h(2),hu(2),hv(2),u(2),v(2),h_hat(2),gamma
    double precision :: alpha(4)
    
    ! Output array intializations
    bmasdq = 0.d0
    bpasdq = 0.d0

    ! Normal and transverse sweep directions
    if (ixy == 1) then
        n_index = 2
        t_index = 3
    else
        n_index = 3
        t_index = 2
    endif
    
    ! ========================================================================
    ! Loop over each cell and decompose fluctuations into transverse waves 
    !  if ixy == 1, and imp == 1 splitting amdq into up and down going
    !  if ixy == 1, and imp == 2 splitting apdq into up and down going
    !  if ixy == 2, and imp == 1 splitting bmdq into left and right going
    !  if ixy == 2, and imp == 2 splitting bpdq into left and right going
    ! ========================================================================
    do i=2-mbc,mx+mbc        
        ! Parse states and pick out important states
        do j=1,2
            layer_index = 3*(j-1)
            ! Solving in the left grid cell (A^-\Delta Q)
            if (imp == 1) then
                h(j) = qr(layer_index + 1,i-1) / rho(j)
                hu(j) = qr(layer_index + n_index,i-1) / rho(j)
                hv(j) = qr(layer_index + t_index,i-1) / rho(j)
                
                b(3) = aux3(1,i-1)
                b(2) = aux2(1,i-1)
                b(1) = aux1(1,i-1)
                
                h_hat(1) = aux2(aux_layer_index,i-1)
                h_hat(2) = aux2(aux_layer_index+1,i-1)
            ! Solving in the right grid cell (A^+ \Delta Q)
            else
                h(j) = ql(layer_index + 1,i) / rho(j)
                hu(j) = ql(layer_index + n_index,i) / rho(j)
                hv(j) = ql(layer_index + t_index,i) / rho(j)
                
                b(3) = aux3(1,i)
                b(2) = aux2(1,i)
                b(1) = aux1(1,i)
                
                h_hat = aux2(aux_layer_index:aux_layer_index+1,i)
            endif
                
            if (h(j) < dry_tolerance(j)) then
                u(j) = 0.d0
                v(j) = 0.d0
            else
                u(j) = hu(j) / h(j)
                v(j) = hv(j) / h(j)
            endif
        enddo
        
        ! ====================================================================
        !  Check for dry states in the bottom layer, 
        !  This is if we are right next to a wall, use single layer solver
        if (h(1) < dry_tolerance(1) .or. h(2) < dry_tolerance(2)) then
            ! Storage for single layer rpt2
            if (h(1) < dry_tolerance(1)) then
                ql_sl = ql(1:3,i) / rho(1)
                qr_sl = qr(1:3,i-1) / rho(1)
                asdq_sl = asdq(1:3,i) / rho(1)
            else if (h(2) < dry_tolerance(2)) then
                ql_sl = ql(4:6,i) / rho(2)
                qr_sl = qr(4:6,i-1) / rho(2)
                asdq_sl = asdq(4:6,i) / rho(2)
            else
                print *, "Invalid dry-state found."
                print *, "  h = ", h
                stop
            end if

            aux1_sl = aux1(1:3,i-1:i)
            aux2_sl = aux2(1:3,i-1:i)
            aux3_sl = aux3(1:3,i-1:i)
            
            ! Call solve
            call rpt2_single_layer(ixy,imp,ql_sl,qr_sl,aux1_sl,aux2_sl, &
                                   aux3_sl,asdq_sl,bmasdq_sl,bpasdq_sl)
            
            if (h(1) < dry_tolerance(1)) then
                bmasdq(1:3,i) = bmasdq_sl * rho(1)
                bpasdq(1:3,i) = bpasdq_sl * rho(1)
            else if (h(2) < dry_tolerance(2)) then
                bmasdq(4:6,i) = bmasdq_sl * rho(2)
                bpasdq(4:6,i) = bpasdq_sl * rho(2)
            end if
            cycle
        endif
        
        ! ====================================================================
        ! Two-layers with no dry states in the vicinity
        ! Compute eigenvector matrix - Linearized system
        if (eigen_method == 1) then
            call eigen_func(h_hat,h_hat,v,v,u,u,t_index,n_index,s,eig_vec)
        else
            call eigen_func(h,h,v,v,u,u,t_index,n_index,s,eig_vec)
        endif

        ! ====================================================================
        !  Solve projection onto eigenvectors - Use LAPACK's dgesv routine
        !    N - (int) - Number of linear equations (6)
        !    NRHS - (int) - Number of right hand sides (1)
        !    A - (dp(6,6)) - Coefficient matrix, in this case eig_vec
        !    LDA - (int) - Leading dimension of A (6)
        !    IPIV - (int(N)) - Pivot indices
        !    B - (dp(LDB,NRHS)) - RHS of equations (delta)
        !    LDB - (int) - Leading dimension of B (6)
        !    INFO - (int) - Status of result
        !  Note that the solution (betas) are in delta after the call
        delta = asdq(:,i) ! This is what we are splitting up
        A = eig_vec ! We need to do this as the return matrix is modified and
                    ! we have to use eig_vec again to compute fwaves
        call dgesv(6,1,A,6,pivot,delta,6,info)
        if (.not.(info == 0)) then
            print "(a,i2,a,i2)","In transverse solver:  ixy=",ixy," imp=",imp
            print "(a,i3)","  Error solving R beta = delta,",info
            print "(a,i3)","  Location: ",i
            print "(a,6d16.8)","  Eigenspeeds: ",s(:)
            print "(a)","  Eigenvectors:"
            do j=1,6
                print "(a,6d16.8)","  ",(eig_vec(j,mw),mw=1,6)
            enddo
            print "(6d16.8)",h,hu,hv
            stop
        endif

        ! Handle lat-long coordinate systems
        if (coordinate_system == 2) then
            if (ixy == 2) then
                dxdcp=(earth_radius*pi/180.d0)
                dxdcm = dxdcp
            else
                if (imp == 1) then
                    dxdcp = earth_radius*pi*cos(aux3(3,i-1))/180.d0
                    dxdcm = earth_radius*pi*cos(aux1(3,i-1))/180.d0
                else
                    dxdcp = earth_radius*pi*cos(aux3(3,i))/180.d0
                    dxdcm = earth_radius*pi*cos(aux1(3,i))/180.d0
                endif
            endif
        else
            dxdcp = 1.d0
            dxdcm = 1.d0
        endif
        
        ! ====================================================================
        !  Determine transverse fluctuations
        !  We also check to see if the fluctuation would enter a dry cell and
        !  skip that split if this is true
        do mw=1,mwaves
            if (s(mw) > 0.d0) then
                if (h(2) + b(2) < b(3)) then
                    bpasdq(1:3,i) = bpasdq(1:3,i) + dxdcp * s(mw) * delta(mw) * eig_vec(1:3,mw)
                else
                    bpasdq(:,i) = bpasdq(:,i) + dxdcp * s(mw) * delta(mw) * eig_vec(:,mw)
                endif
            else if (s(mw) < 0.d0) then
                if (h(2) + b(2) < b(1)) then
                    bmasdq(1:3,i) = bmasdq(1:3,i) + dxdcm * s(mw) * delta(mw) * eig_vec(1:3,mw)
                else
                    bmasdq(:,i) = bmasdq(:,i) + dxdcm * s(mw) * delta(mw) * eig_vec(:,mw)
                endif
            endif
        enddo
    enddo

end subroutine rpt2


subroutine rpt2_single_layer(ixy,imp,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! Single layer point-wise transverse Riemann solver using an einfeldt Jacobian
! Note that there have been some changes to variable definitions in this
! routine from the original vectorized one. 
!
! Adapted from geoclaw 4-23-2011

    use amr_module, only: mcapa

    use geoclaw_module, only: g => grav, earth_radius, pi, coordinate_system
    use multilayer_module, only: num_layers, eigen_method, inundation_method
    use multilayer_module, only: dry_tolerance

    implicit none

    integer, parameter :: meqn = 3
    integer, parameter :: mwaves = 3

    integer, intent(in) :: ixy,imp

    real(kind=8), intent(in) :: ql(meqn) ! = ql(i,meqn)
    real(kind=8), intent(in) :: qr(meqn) ! = qr(i-1,meqn)
    real(kind=8), intent(in out) :: asdq(meqn)
    real(kind=8), intent(in out) :: bmasdq(meqn)
    real(kind=8), intent(in out) :: bpasdq(meqn)
    ! Since we need two values of aux, 1 = i-1 and 2 = i
    real(kind=8), intent(in) :: aux1(3,2)
    real(kind=8), intent(in) :: aux2(3,2)
    real(kind=8), intent(in) :: aux3(3,2)

    real(kind=8) :: s(3)
    real(kind=8) :: r(3,3)
    real(kind=8) :: beta(3)
    real(kind=8) :: abs_tol
    real(kind=8) :: hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
    real(kind=8) :: uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
    real(kind=8) :: delf1,delf2,delf3,dxdcd,dxdcu
    real(kind=8) :: dxdcm,dxdcp,topo1,topo3,eta,tol

    integer ::  m,mw,mu,mv
    
    tol = dry_tolerance(1)
    abs_tol = tol

    if (ixy.eq.1) then
     mu = 2
     mv = 3
    else
     mu = 3
     mv = 2
    endif


       hl=qr(1)
       hr=ql(1)
       hul=qr(mu)
       hur=ql(mu)
       hvl=qr(mv)
       hvr=ql(mv)

!===========determine velocity from momentum===========================
       if (hl.lt.abs_tol) then
          hl=0.d0
          ul=0.d0
          vl=0.d0
       else
          ul=hul/hl
          vl=hvl/hl
       endif

       if (hr.lt.abs_tol) then
          hr=0.d0
          ur=0.d0
          vr=0.d0
       else
          ur=hur/hr
          vr=hvr/hr
       endif

       do mw=1,mwaves
          s(mw)=0.d0
          beta(mw)=0.d0
          do m=1,meqn
             r(m,mw)=0.d0
          enddo
       enddo
      dxdcp = 1.d0
      dxdcm = 1.d0

       if (hl.le.dry_tolerance(1).and.hr.le.dry_tolerance(1)) go to 90

       !check and see if cell that transverse waves are going in is high and dry
       if (imp.eq.1) then
            eta = qr(1) + aux2(1,1)
            topo1 = aux1(1,1)
            topo3 = aux3(1,1)
       else
            eta = ql(1) + aux2(1,2)
            topo1 = aux1(1,2)
            topo3 = aux3(1,2)
       endif
       if (eta.lt.max(topo1,topo3)) go to 90

      if (coordinate_system == 2) then
         if (ixy.eq.2) then
            dxdcp=(earth_radius*pi/180.d0)
            dxdcm = dxdcp
         else
            if (imp.eq.1) then
               dxdcp = earth_radius*pi*cos(aux3(3,1))/180.d0
               dxdcm = earth_radius*pi*cos(aux1(3,1))/180.d0
            else
               dxdcp = earth_radius*pi*cos(aux3(3,2))/180.d0
               dxdcm = earth_radius*pi*cos(aux1(3,2))/180.d0
            endif
         endif
      endif


!=====Determine some speeds necessary for the Jacobian=================
            vhat=(vr*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + &
             (vl*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))

            uhat=(ur*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + &
             (ul*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
            hhat=(hr+hl)/2.d0

            roe1=vhat-dsqrt(g*hhat)
            roe3=vhat+dsqrt(g*hhat)

            s1l=vl-dsqrt(g*hl)
            s3r=vr+dsqrt(g*hr)

            s1=dmin1(roe1,s1l)
            s3=dmax1(roe3,s3r)

            s2=0.5d0*(s1+s3)

           s(1)=s1
           s(2)=s2
           s(3)=s3
!=======================Determine asdq decomposition (beta)============
         delf1=asdq(1)
         delf2=asdq(mu)
         delf3=asdq(mv)

         beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
         beta(2) = -s2*delf1 + delf2
         beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
!======================End =================================================

!=====================Set-up eigenvectors===================================
         r(1,1) = 1.d0
         r(2,1) = s2
         r(3,1) = s1

         r(1,2) = 0.d0
         r(2,2) = 1.d0
         r(3,2) = 0.d0

         r(1,3) = 1.d0
         r(2,3) = s2
         r(3,3) = s3
!============================================================================
90      continue
!============= compute fluctuations==========================================

            do  m=1,meqn
               bmasdq(m)=0.0d0
               bpasdq(m)=0.0d0
            enddo
            do  mw=1,3
               if (s(mw).lt.0.d0) then
                 bmasdq(1) =bmasdq(1) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                 bmasdq(mu)=bmasdq(mu)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                 bmasdq(mv)=bmasdq(mv)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
               elseif (s(mw).gt.0.d0) then
                 bpasdq(1) =bpasdq(1) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                 bpasdq(mu)=bpasdq(mu)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                 bpasdq(mv)=bpasdq(mv)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
               endif
            enddo


end subroutine rpt2_single_layer

