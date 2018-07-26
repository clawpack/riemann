module acoustics_module
    contains
    ! ==============================================================================
    !  Point-wise version of the homogeneous 2d acoustics problem
    ! ==============================================================================

    attributes(device) &
    subroutine riemann_acoustics_homo_2d(ixy, num_eqn, num_aux, num_waves, q_l, q_r, aux_l, aux_r,    &
                          wave, s, bulk, rho) &
                          bind (C,name='riemann_acoustics_homo_2d_fortran')

        implicit none

        ! Input Arguments
        integer, value, intent(in) :: ixy, num_eqn, num_aux, num_waves
        real(CLAW_REAL), value, intent(in) :: rho, bulk
        real(CLAW_REAL), intent(in) :: q_l(num_eqn), q_r(num_eqn)
        real(CLAW_REAL), intent(in) :: aux_l(num_aux), aux_r(num_aux)

        ! Output arguments
        real(CLAW_REAL), intent(out) :: wave(num_eqn, num_waves)
        real(CLAW_REAL), intent(out) :: s(num_waves)

        ! Locals
        integer :: normal_index, transverse_index
        real(CLAW_REAL) :: delta(3), a(2)
        real(CLAW_REAL) :: c,z


        c = sqrt(bulk/rho)
        z = c*rho

        if (ixy == 1) then
            normal_index = 2
            transverse_index = 3
        else
            normal_index = 3
            transverse_index = 2
        end if

        ! note that notation for u and v reflects assumption that the
        ! Riemann problems are in the x-direction with u in the normal
        ! direciton and v in the orthogonal direcion, but with the above
        ! definitions of mu and mv the routine also works with ixy=2
        ! in which case waves come from the
        ! Riemann problems u_t + g(u)_y = 0 in the y-direction.


        ! Split the jump in q at each interface into waves
        delta = q_r - q_l
        ! Find coefficients for the 2 eigenvectors:
        a(1) = (-delta(1) + z * delta(normal_index)) / (2.d0 * z)
        a(2) = ( delta(1) + z * delta(normal_index)) / (2.d0 * z)

        ! Compute waves
        wave(1, 1) = -a(1) * z
        wave(normal_index, 1) = a(1)
        wave(transverse_index, 1) = 0.d0
        s(1) = -c

        wave(1, 2) = a(2) * z
        wave(normal_index, 2) = a(2)
        wave(transverse_index, 2) = 0.d0
        s(2) = c

    end subroutine riemann_acoustics_homo_2d

end module acoustics_module
