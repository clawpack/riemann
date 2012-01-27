

!     =====================================================
    subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx, &
                ql,qr,aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq,num_aux)
!     =====================================================
    implicit double precision (a-h,o-z)

!     # Dummy transverse Riemann solver, for use in dimensionally-split algorithm.

    write(*,*) 'Error: Dummy transverse Riemann solver called!'
    return
    end subroutine rpt2
