#include "arrayIndex.H"

#define GRAVITY 1.d0

module shallow_module
    contains
    ! =====================================================
    attributes(device) &
    subroutine riemann_shallow_x(q_l, q_r, auxl, auxr,    &
                          wave, s) &
                          bind (C,name='riemann_shallow_x')
    ! =====================================================

        ! Roe-solver for the 2D shallow water equations
        ! solve Riemann problems along one slice of data.

        ! waves: 3
        ! equations: 3

        ! Conserved quantities:
        !       1 depth
        !       2 x_momentum
        !       3 y_momentum

        implicit none

        real(CLAW_REAL), intent(in) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: auxl(NCOEFFS), auxr(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: wave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        !   locals
        !   ------------
        real(CLAW_REAL) :: u,v,a,h
        real(CLAW_REAL) :: hsqrtl, hsqrtr, hsq2
        real(CLAW_REAL) :: a1,a2,a3
        real(CLAW_REAL) :: delta(3)

        !   # Compute the Roe-averaged variables needed in the Roe solver.

        h = (q_l(1)+q_r(1))*0.50d0
        hsqrtl = dsqrt(q_l(1))
        hsqrtr = dsqrt(q_r(1))
        hsq2 = hsqrtl + hsqrtr
        u = (q_l(2)/hsqrtl + q_r(2)/hsqrtr) / hsq2
        v = (q_l(3)/hsqrtl + q_r(3)/hsqrtr) / hsq2
        a = dsqrt(GRAVITY*h)


        !   # now split the jump in q at each interface into waves
        !   # find a1 thru a3, the coefficients of the 3 eigenvectors:
        delta(1) = q_r(1) - q_l(1)
        delta(2) = q_r(2) - q_l(2)
        delta(3) = q_r(3) - q_l(3)
        a1 = ((u+a)*delta(1) - delta(2))*(0.50d0/a)
        a2 = -v*delta(1) + delta(3)
        a3 = (-(u-a)*delta(1) + delta(2))*(0.50d0/a)

        !      # Compute the waves.

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1*(u-a)
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1*v
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 1, NWAVES, blockDim%y, blockDim%x) &
        ) = u-a

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = 0.0d0
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = 0.0d0
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a2
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 2, NWAVES, blockDim%y, blockDim%x) &
        ) = u

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3*(u+a)
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3*v
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 3, NWAVES, blockDim%y, blockDim%x) &
        ) = u+a
        return
    end subroutine riemann_shallow_x


    attributes(device) &
    subroutine riemann_shallow_y(q_l, q_r, auxl, auxr,    &
                          wave, s) &
                          bind (C,name='riemann_shallow_y')
    ! =====================================================

        ! Roe-solver for the 2D shallow water equations
        ! solve Riemann problems along one slice of data.

        ! waves: 3
        ! equations: 3

        ! Conserved quantities:
        !       1 depth
        !       2 x_momentum
        !       3 y_momentum

        implicit none

        real(CLAW_REAL), intent(in) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: auxl(NCOEFFS), auxr(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: wave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        !   locals
        !   ------------
        real(CLAW_REAL) :: u,v,a,h
        real(CLAW_REAL) :: hsqrtl, hsqrtr, hsq2
        real(CLAW_REAL) :: a1,a2,a3
        real(CLAW_REAL) :: delta(3)

        !   # Compute the Roe-averaged variables needed in the Roe solver.

        h = (q_l(1)+q_r(1))*0.50d0
        hsqrtl = dsqrt(q_l(1))
        hsqrtr = dsqrt(q_r(1))
        hsq2 = hsqrtl + hsqrtr
        u = (q_l(3)/hsqrtl + q_r(3)/hsqrtr) / hsq2
        v = (q_l(2)/hsqrtl + q_r(2)/hsqrtr) / hsq2
        a = dsqrt(GRAVITY*h)


        !   # now split the jump in q at each interface into waves
        !   # find a1 thru a3, the coefficients of the 3 eigenvectors:
        delta(1) = q_r(1) - q_l(1)
        delta(2) = q_r(3) - q_l(3)
        delta(3) = q_r(2) - q_l(2)
        a1 = ((u+a)*delta(1) - delta(2))*(0.50d0/a)
        a2 = -v*delta(1) + delta(3)
        a3 = (-(u-a)*delta(1) + delta(2))*(0.50d0/a)

        !      # Compute the waves.

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1*(u-a)
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a1*v
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 1, NWAVES, blockDim%y, blockDim%x) &
        ) = u-a

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = 0.0d0
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = 0.0d0
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a2
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 2, NWAVES, blockDim%y, blockDim%x) &
        ) = u

        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3*(u+a)
        wave(&
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 3, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
        ) = a3*v
        s(&
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 3, NWAVES, blockDim%y, blockDim%x) &
        ) = u+a
        return
    end subroutine riemann_shallow_y
end module shallow_module

! This is for the example, GPU/shallow_water_no_topo
! These subroutines are not acutally used, but only for compilation
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
end subroutine rpn2

subroutine rpt2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
end subroutine rpt2
