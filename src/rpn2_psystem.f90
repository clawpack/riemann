! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================

! Approximate Riemann solver for the nonlinear p-system in 2d
! with variable coefficients.

! waves: 2
! equations: 3

! Conserved quantities:
!       1 strain
!       2 x_momentum
!       3 y_momentum

! The variable names are based on an elasticity interpretation

! The jacobian matrix of the flux vector (in each direction)
! is approximated by the linear localized problem.
 
! There are 3 eigenvectors; however, the second eigenvalue
! is always zero and just two waves are computed.

! Solve Riemann problems along one slice of data:
!  in the x-direction if ixy=1
!  in the y-direction if ixy=2.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! f-wave approach is considered. This consists in decomposing the flux
! difference (assuming it's continuous) using the eigenvectors of the
! jacobian matrices of the flux vectors (in each spatial direction).
! This is convinient for variable coefficient PDEs.

! On output, fwave contains the f-waves,
!            s the speeds,
!            amdq the  left-going fluctuation
!            apdq the right-going fluctuation

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit double precision (a-h,o-z)
    dimension fwave(meqn,mwaves,1-mbc:maxm+mbc)
    dimension    s(mwaves,1-mbc:maxm+mbc)
    dimension   ql(meqn,1-mbc:maxm+mbc)
    dimension   qr(meqn,1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)
    dimension auxl(maux,1-mbc:maxm+mbc)
    dimension auxr(maux,1-mbc:maxm+mbc)
          
    do 10 i=2-mbc,mx+mbc
    ! material properties
        pim=auxr(1,i-1)
        pi=auxl(1,i)
        Eim=auxr(2,i-1)
        Ei=auxl(2,i)
    ! solution eps, urho and vrho
        epsi=ql(1,i)
        urhoi=ql(2,i)
        vrhoi=ql(3,i)
        epsim=qr(1,i-1)
        urhoim=qr(2,i-1)
        vrhoim=qr(3,i-1)
    ! linearity of material (for cell i and for cell im)
        linearity_mati=auxl(3,i)
        linearity_matim=auxr(3,i-1)
    ! sigma
        sigmai=sigma(epsi,Ei,linearity_mati)
        sigmaim=sigma(epsim,Eim,linearity_matim)
    ! sigmap
        sigmapi=sigmap(epsi,Ei,linearity_mati)
        sigmapim=sigmap(epsim,Eim,linearity_matim)
    ! computation of components of eigenvectors
        r11=1/dsqrt(sigmapim*pim)
        r13=-1/dsqrt(sigmapi*pi)
    ! shock speeds
        s(1,i)=-dsqrt(sigmapim/pim)  !lambda_1
        s(2,i)=dsqrt(sigmapi/pi)     !lambda_2
             
        if(ixy == 1) then      !x direction
        ! compute jump in flux
            dF1=-(urhoi/pi-urhoim/pim)
            dF2=-(sigmai-sigmaim)
        ! compute betas
            beta1=(dF1-r13*dF2)/(r11-r13)
            beta3=(-dF1+r11*dF2)/(r11-r13)
        ! compute f-waves
            fwave(1,1,i)=beta1*r11
            fwave(2,1,i)=beta1*1
            fwave(3,1,i)=beta1*0
            fwave(1,2,i)=beta3*r13
            fwave(2,2,i)=beta3*1
            fwave(3,2,i)=beta3*0
        else                   !y direction
        ! compute jump in flux
            dF1=-(vrhoi/pi-vrhoim/pim)
            dF3=-(sigmai-sigmaim)
        ! compute betas
            beta1=(dF1-r13*dF3)/(r11-r13)
            beta3=(-dF1+r11*dF3)/(r11-r13)
        ! compute f-waves
            fwave(1,1,i)=beta1*r11
            fwave(2,1,i)=beta1*0
            fwave(3,1,i)=beta1*1
            fwave(1,2,i)=beta3*r13
            fwave(2,2,i)=beta3*0
            fwave(3,2,i)=beta3*1
        endif
    ! computation of the fluctuations
        amdq(1,i)=fwave(1,1,i)
        amdq(2,i)=fwave(2,1,i)
        amdq(3,i)=fwave(3,1,i)
        apdq(1,i)=fwave(1,2,i)
        apdq(2,i)=fwave(2,2,i)
        apdq(3,i)=fwave(3,2,i)

    10 END DO
         
    return
    end subroutine rpn2

    double precision function sigma(eps,E,linearity_mat)
!     Returns the flux sigma for a given
!     eps, E and depending the linearity of the material
    implicit double precision (a-h,o-z)
    beta=5
    select case (linearity_mat)
    case (1)
    sigma=E*eps
    case (2)
    sigma=dexp(E*eps)-1
    case (3)
    sigma=0.1*E*eps+beta*eps**3*E**3
    end select
    return
    END function

    double precision function sigmap(eps,E,linearity_mat)
    implicit double precision (a-h,o-z)
! Returns the derivative of sigma wrt eps for a given
!     eps, E and depending the linearity of the material
    beta=5
    select case (linearity_mat)
    case (1)
    sigmap=E
    case (2)
    sigmap=E*dexp(E*eps)
    case (3)
    sigmap=0.1*E+3*beta*eps**2*E**3
    end select
    return
    END function
