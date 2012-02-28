module rp1_layered_shallow_water_module

    implicit none
    save

    ! Precision    
    integer, private, parameter :: R = kind(1.d0)
    
    ! Method parameters
    integer :: eigen_method,inundation_method
    real(kind=R) :: dry_tolerance
    logical :: entropy_fix

end module rp1_layered_shallow_water_module