!-----------------------------------------------------------------------------------
! Name: rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,qwave,s,amdq,apdq)
!
! Description: Riemann Q-wave solver for the 3D Euler Equations
!
! Method: Riemann solver for 3D Euler (neutral fluid) with Q-Wave jumps.
!         Solves along 1D slices of data (uses dimensional splitting).
!         Roe-averaged quantities are used for the 'cell centered'
!         quantities.
!
! Inputs: 
!         ixyz <INTEGER>   : direction to take slice x-direction if ixyz=1
!                                                    y-direction if ixyz=2.
!                                                    z-direction if ixyz=3.
!         maxm <INTEGER>   : max number of grid cells (less the ghost cells)
!         meqn <INTEGER>   : number of equations in the system
!         mwaves <INTEGER> : nuber of waves in the system
!         maux <INTEGER>   : number of auxilary equations
!         mbc <INTEGER>    : number of ghost cells on either end
!         mx <INTEGER>     : number of elements
!         ql <REAL>        : state vector at left edge of each cell
!                            Note the i'th Riemann problem has left state qr(:,i-1)
!         qr <REAL>        : state vector at right edge of each cell
!                            Note the i'th Riemann problem has right state ql(:,i)
!         auxl <REAL>      : state of auxilary variable on left egde of cell 
!         auxr <REAL>      : state of auxilary variable on right egde of cell 
!         
! Outputs: 
!          qwave <REAL>    : q-wave vectors
!          s <REAL>        : eigenvalues (wave speeds)
!          amdq <REAL>     : left-going flux differences
!          apdq <REAL>     : right-going flux differences
!
!-----------------------------------------------------------------------------------
SUBROUTINE rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,qwave,s,amdq,apdq)

  !-----------------------------------------------------------------------------------
  ! Declarations
  !-----------------------------------------------------------------------------------
  IMPLICIT NONE

  ! Input
  INTEGER, INTENT(IN) :: ixyz,maxm,meqn,mwaves,maux,mbc,mx
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(IN) :: ql,qr,auxl,auxr

  ! Output
  REAL(kind=8), DIMENSION(meqn,mwaves,1-mbc:maxm+mbc), INTENT(INOUT) :: qwave
  REAL(kind=8), DIMENSION(mwaves,1-mbc:maxm+mbc), INTENT(INOUT) :: s
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(INOUT) :: amdq,apdq
  
  ! Local Storage
  INTEGER, PARAMETER :: maxmrp = 1002
  REAL(kind=8), DIMENSION(5) :: delta,alpha
  REAL(kind=8), DIMENSION(-1:maxmrp) :: u2v2w2,u,v,w,enth,a,g1a2,euv
  REAL(kind=8) :: asqrd,c1,c2,ci,cim1,en1,en2,p1,p2,pi,pim1,pl,pr,&
       rho1,rho2,rhoi,rhoim1,rhou1,rhou2,rhov1,rhov2,rhow1,rhow2,&
       rhsq2,rhsqrtl,rhsqrtr,s0,s1,s2,s3,sfract,df
  INTEGER :: i,j,m,mu,mv,mw,mws

  ! Use entropy fix for transonic rarefactions
  LOGICAL, PARAMETER :: efix = .TRUE.
  
  ! Common block storage for ideal gas constant
  REAL(kind=8) :: gamma,gamma1
  COMMON /cparam/ gamma

  ! Set (gamma-1)
  gamma1 = gamma - 1.d0


  !-----------------------------------------------------------------------------------
  ! Make sure maxmrp is large enough
  !-----------------------------------------------------------------------------------
  !if (-1 > 1-mbc .or. maxmrp < maxm+mbc) then
  IF (maxmrp < maxm+mbc)THEN
     WRITE(6,*) 'need to increase maxmrp in rpn3_neutral_qwave.f90'
     WRITE(6,*) 'maxmrp: ',maxmrp,' maxm: ',maxm,' mbc: ',mbc
     WRITE(6,*) 'maxm+mbc=',maxm+mbc
     STOP
  ENDIF

  !-----------------------------------------------------------------------------------
  ! Make sure the number of waves is 3 for the 3D Euler equations
  !-----------------------------------------------------------------------------------
  IF (mwaves /= 3) THEN
     WRITE(6,*) '*** Should have mwaves=3 for this Riemann solver'
     STOP
  ENDIF
    
  !-----------------------------------------------------------------------------------
  ! Set mu to point to  the component of the system that corresponds to momentum in 
  !  the direction of this slice, mv and mw to the orthogonal momentums:
  !-----------------------------------------------------------------------------------
  IF(ixyz == 1)THEN
     mu = 2
     mv = 3
     mw = 4
  ELSE IF(ixyz == 2)THEN
     mu = 3
     mv = 4
     mw = 2
  ELSE
     mu = 4
     mv = 2
     mw = 3
  ENDIF

  !---------------------------------------------------------------------------------
  ! Note that notation for u,v, and w reflects assumption that the
  !   Riemann problems are in the x-direction with u in the normal
  !   direction and v and w in the orthogonal directions, but with the
  !   above definitions of mu, mv, and mw the routine also works with
  !   ixyz=2 and ixyz = 3
  !   and returns, for example, f0 as the Godunov flux g0 for the
  !   Riemann problems u_t + g(u)_y = 0 in the y-direction.
  !
  ! Compute the Roe-averaged variables needed in the Roe solver.
  !---------------------------------------------------------------------------------
  
  !---------------------------------------------------------------------------------
  ! Initialize qwaves and speeds
  !---------------------------------------------------------------------------------
  qwave = 0.d0
  s = 0.d0

  !---------------------------------------------------------------------------------
  ! Loop over grid cell interfaces
  !---------------------------------------------------------------------------------
  DO i = 2-mbc, mx+mbc
     IF (qr(1,i-1) <= 0.d0 .OR. ql(1,i) <= 0.d0) THEN
        WRITE(*,*) i, mu, mv, mw
        WRITE(*,'(5e12.4)') (qr(j,i-1),j=1,5)
        WRITE(*,'(5e12.4)') (ql(j,i),j=1,5)
        IF (ixyz == 1) &
             WRITE(6,*) '*** rho <= 0 in x-sweep at ',i
        IF (ixyz == 2) &
             WRITE(6,*) '*** rho <= 0 in y-sweep at ',i
        IF (ixyz == 3) &
             WRITE(6,*) '*** rho <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with rho <= 0...'
        STOP
     ENDIF
     rhsqrtl = SQRT(qr(1,i-1))
     rhsqrtr = SQRT(ql(1,i))
     pl = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 + &
          qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1))
     pr = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 + &
          ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i))
     rhsq2 = rhsqrtl + rhsqrtr
     u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
     v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
     w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
     enth(i) = (((qr(5,i-1)+pl)/rhsqrtl &
          + (ql(5,i)+pr)/rhsqrtr)) / rhsq2
     u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
     asqrd = gamma1*(enth(i) - .5d0*u2v2w2(i))
     
     IF (i>=0 .AND. i<=mx .AND. asqrd <= 0.d0) THEN
        IF (ixyz == 1) &
             WRITE(6,*) '*** a**2 <= 0 in x-sweep at ',i
        IF (ixyz == 2) &
             WRITE(6,*) '*** a**2 <= 0 in y-sweep at ',i
        IF (ixyz == 3) &
             WRITE(6,*) '*** a**2 <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with a**2 < 0...'
        STOP
     ENDIF
     a(i) = SQRT(asqrd)
     g1a2(i) = gamma1 / asqrd
     euv(i) = enth(i) - u2v2w2(i)
  END DO

  !---------------------------------------------------------------------------------
  ! Now split the jump in q1d at each interface into q-waves
  !   find alpha(1) thru alpha(5), the coefficients of the 5 eigenvectors:
  !---------------------------------------------------------------------------------
  DO i = 2-mbc, mx+mbc
     delta(1) = ql(1,i) - qr(1,i-1)
     delta(2) = ql(mu,i) - qr(mu,i-1)
     delta(3) = ql(mv,i) - qr(mv,i-1)
     delta(4) = ql(mw,i) - qr(mw,i-1)
     delta(5) = ql(5,i) - qr(5,i-1)
     alpha(4) = g1a2(i) * (euv(i)*delta(1) &
          + u(i)*delta(2) + v(i)*delta(3) + w(i)*delta(4) - delta(5))
     alpha(2) = delta(3) - v(i)*delta(1)
     alpha(3) = delta(4) - w(i)*delta(1)
     alpha(5) = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*alpha(4)) / (2.d0*a(i))
     alpha(1) = delta(1) - alpha(4) - alpha(5)
     
     !------------------------------------------------------------------------------
     ! Compute the q-waves.
     !   Note that the 2-wave, 3-wave and 4-wave travel at the same speed
     !   and are lumped together in qwave(.,2,.).  The 5-wave is then stored
     !   in qwave(.,3,.).
     !------------------------------------------------------------------------------
     qwave(1,1,i)  = alpha(1)
     qwave(mu,1,i) = alpha(1)*(u(i)-a(i))
     qwave(mv,1,i) = alpha(1)*v(i)
     qwave(mw,1,i) = alpha(1)*w(i)
     qwave(5,1,i)  = alpha(1)*(enth(i) - u(i)*a(i))
     s(1,i) = u(i)-a(i)
     
     qwave(1,2,i)  = alpha(4)
     qwave(mu,2,i) = alpha(4)*u(i)
     qwave(mv,2,i) = alpha(4)*v(i) + alpha(2)
     qwave(mw,2,i) = alpha(4)*w(i) + alpha(3)
     qwave(5,2,i)  = alpha(4)*0.5d0*u2v2w2(i)  + alpha(2)*v(i) + alpha(3)*w(i)
     s(2,i) = u(i)
     
     qwave(1,3,i)  = alpha(5)
     qwave(mu,3,i) = alpha(5)*(u(i)+a(i))
     qwave(mv,3,i) = alpha(5)*v(i)
     qwave(mw,3,i) = alpha(5)*w(i)
     qwave(5,3,i)  = alpha(5)*(enth(i)+u(i)*a(i))
     s(3,i) = u(i)+a(i)
  END DO

  !---------------------------------------------------------------------------------
  !    # compute flux differences amdq and apdq.
  !---------------------------------------------------------------------------------
  IF (.NOT.efix) THEN
     !------------------------------------------------------------------------------
     !  No entropy fix
     !------------------------------------------------------------------------------
     !  amdq = SUM s*qwave   over left-going waves
     !  apdq = SUM s*qwave   over right-going waves
     !------------------------------------------------------------------------------
     amdq = 0.d0
     apdq = 0.d0
     DO i=2-mbc, mx+mbc
        DO mws=1,mwaves
           IF (s(mws,i) < 0.d0) THEN
              amdq(:,i) = amdq(:,i) + s(mws,i)*qwave(:,mws,i)
           ELSE
              apdq(:,i) = apdq(:,i) + s(mws,i)*qwave(:,mws,i)
           ENDIF
        END DO
     END DO

  ELSE
     !------------------------------------------------------------------------------
     ! With entropy fix
     !------------------------------------------------------------------------------
     !------------------------------------------------------------------------------
     ! Compute flux differences amdq and apdq.
     !   First compute amdq as sum of s*qwave for left going waves.
     !   Incorporate entropy fix by adding a modified fraction of wave
     !   if s should change sign.
     !------------------------------------------------------------------------------
     amdq = 0.d0
     apdq = 0.d0
     DO i = 2-mbc, mx+mbc

        !---------------------------------------------------------------------------
        ! Check 1-wave
        !---------------------------------------------------------------------------
        rhoim1 = qr(1,i-1)
        pim1 = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 &
             + qr(mv,i-1)**2 + qr(mw,i-1)**2) / rhoim1)
        cim1 = SQRT(gamma*pim1/rhoim1)
        s0 = qr(mu,i-1)/rhoim1 - cim1  ! u-c in left state (cell i-1)

        ! Check for fully supersonic case:
        IF (s0 >= 0.d0 .AND. s(1,i) > 0.d0)CYCLE ! everything is right-going

        rho1 = qr(1,i-1) + qwave(1,1,i)
        rhou1 = qr(mu,i-1) + qwave(mu,1,i)
        rhov1 = qr(mv,i-1) + qwave(mv,1,i)
        rhow1 = qr(mw,i-1) + qwave(mw,1,i)
        en1 = qr(5,i-1) + qwave(5,1,i)
        p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 + rhow1**2)/rho1)
        c1 = SQRT(gamma*p1/rho1)
        s1 = rhou1/rho1 - c1  ! u-c to right of 1-wave
        IF (s0 < 0.d0 .AND. s1 > 0.d0) THEN
           ! transonic rarefaction in the 1-wave
           sfract = s0 * (s1-s(1,i)) / (s1-s0)
        ELSE IF (s(1,i) < 0.d0) THEN
           ! 1-wave is leftgoing
           sfract = s(1,i)
        ELSE
           ! 1-wave is rightgoing
           sfract = 0.d0   ! this shouldn't happen since s0 < 0
        ENDIF
        amdq(:,i) = sfract*qwave(:,1,i)
           
        !---------------------------------------------------------------------------
        ! Check 2-wave
        !---------------------------------------------------------------------------
        IF (s(2,i) >= 0.d0)CYCLE ! 2-,3- and 4- waves are rightgoing
        amdq(:,i) = amdq(:,i) + s(2,i)*qwave(:,2,i)
        
        !---------------------------------------------------------------------------
        ! Check 3-wave
        !---------------------------------------------------------------------------
        rhoi = ql(1,i)
        pi = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2+ql(mv,i)**2+ql(mw,i)**2)/rhoi)
        ci = SQRT(gamma*pi/rhoi)
        s3 = ql(mu,i)/rhoi + ci     ! u+c in right state  (cell i)
        
        rho2 = ql(1,i) - qwave(1,3,i)
        rhou2 = ql(mu,i) - qwave(mu,3,i)
        rhov2 = ql(mv,i) - qwave(mv,3,i)
        rhow2 = ql(mw,i) - qwave(mw,3,i)
        en2 = ql(5,i) - qwave(5,3,i)
        p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 + rhow2**2)/rho2)
        c2 = SQRT(gamma*p2/rho2)
        s2 = rhou2/rho2 + c2   !# u+c to left of 3-wave
        IF (s2 < 0.d0 .AND. s3 > 0.d0 ) THEN
           ! transonic rarefaction in the 3-wave
           sfract = s2 * (s3-s(3,i)) / (s3-s2)
        ELSE IF (s(3,i) < 0.d0) THEN
           ! 3-wave is leftgoing
           sfract = s(3,i)
        ELSE
           ! 3-wave is rightgoing
           CYCLE
        ENDIF
        
        amdq(:,i) = amdq(:,i) + sfract*qwave(:,3,i)
     END DO
     
     !---------------------------------------------------------------------------
     ! Compute the remaining right-going flux differences:
     ! df = SUM s*wave   is the total flux difference and apdq = df - amdq
     !---------------------------------------------------------------------------
     DO m=1,meqn
        DO i = 2-mbc, mx+mbc
           df = 0.d0
           DO mws=1,mwaves
              df = df + s(mws,i)*qwave(m,mws,i)
           END DO
           apdq(m,i) = df - amdq(m,i)
        END DO
     END DO

  END IF ! efix loop

  !--------------------------------------------------------------------------------
  ! Terminate
  !--------------------------------------------------------------------------------
END SUBROUTINE rpn3
