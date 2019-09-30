#include "arrayIndex.H"

module acoustics_module
    use setprob_module, only: bulk => bulk_d, rho => rho_d
    implicit none
    contains
    ! ==============================================================================
    !  Point-wise version of the homogeneous 2d acoustics problem
    ! ==============================================================================

    attributes(device) &
    subroutine riemann_acoustics_homo_2d_x(q_l, q_r, aux_l, aux_r,    &
                          wave, s) &
                          bind (C,name='riemann_acoustics_homo_2d_fortran_x')

        implicit none

        real(CLAW_REAL), intent(in) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: aux_l(NCOEFFS), aux_r(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: wave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        ! Locals
        real(CLAW_REAL) :: delta(3), a(2)
        real(CLAW_REAL) :: c,z


        c = sqrt(bulk/rho)
        z = c*rho

        ! note that notation for u and v reflects assumption that the
        ! Riemann problems are in the x-direction with u in the normal
        ! direciton and v in the orthogonal direcion, but with the above
        ! definitions of mu and mv the routine also works with ixy=2
        ! in which case waves come from the
        ! Riemann problems u_t + g(u)_y = 0 in the y-direction.


        ! Split the jump in q at each interface into waves
        delta = q_r - q_l
        ! Find coefficients for the 2 eigenvectors:
        a(1) = (-delta(1) + z * delta(2)) / (2.d0 * z)
        a(2) = ( delta(1) + z * delta(2)) / (2.d0 * z)

        ! Compute waves
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = -a(1) * z
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(1)
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = 0.d0
        s( &
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 1, NWAVES, blockDim%y, blockDim%x) &
            ) = -c

        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(2) * z
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(2)
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = 0.d0
        s( &
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 2, NWAVES, blockDim%y, blockDim%x) &
            ) = c


    end subroutine riemann_acoustics_homo_2d_x

    attributes(device) &
    subroutine riemann_acoustics_homo_2d_y(q_l, q_r, aux_l, aux_r,    &
                          wave, s) &
                          bind (C,name='riemann_acoustics_homo_2d_fortran_y')

        implicit none

        real(CLAW_REAL), intent(in) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: aux_l(NCOEFFS), aux_r(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: wave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        ! Locals
        real(CLAW_REAL) :: delta(3), a(2)
        real(CLAW_REAL) :: c,z


        c = sqrt(bulk/rho)
        z = c*rho

        ! note that notation for u and v reflects assumption that the
        ! Riemann problems are in the x-direction with u in the normal
        ! direciton and v in the orthogonal direcion, but with the above
        ! definitions of mu and mv the routine also works with ixy=2
        ! in which case waves come from the
        ! Riemann problems u_t + g(u)_y = 0 in the y-direction.


        ! Split the jump in q at each interface into waves
        delta = q_r - q_l
        ! Find coefficients for the 2 eigenvectors:
        a(1) = (-delta(1) + z * delta(3)) / (2.d0 * z)
        a(2) = ( delta(1) + z * delta(3)) / (2.d0 * z)

        ! Compute waves

        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = -a(1) * z
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(1)
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 1, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = 0.d0
        s( &
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 1, NWAVES, blockDim%y, blockDim%x) &
            ) = -c

        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(2) * z
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = a(2)
        wave( &
            GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, 2, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
            ) = 0.d0
        s( &
            GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, 2, NWAVES, blockDim%y, blockDim%x) &
            ) = c



    end subroutine riemann_acoustics_homo_2d_y
end module acoustics_module
