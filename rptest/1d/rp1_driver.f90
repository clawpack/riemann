
module rp1_driver

contains

    subroutine call_rp1(meqn,mwaves,q_left,q_right,wave,s,amdq,apdq)
    ! wrapper to call a 1d Riemann solver
    ! still need to add auxl and auxr

    implicit none

    ! Input:
    integer, intent(in) :: meqn, mwaves
    real(kind=8), intent(in) ::   q_left(meqn)
    real(kind=8), intent(in) ::   q_right(meqn)

    ! Output:
    real(kind=8), intent(out) ::    s(mwaves)
    real(kind=8), intent(out) :: wave(meqn, mwaves)
    real(kind=8), intent(out) ::  apdq(meqn)
    real(kind=8), intent(out) ::  amdq(meqn)

    ! local variables:
    integer, parameter :: mbc = 0
    integer, parameter :: maux = 0    ! NEED TO ADD aux arrays
    real(kind=8) ::   ql(meqn, 1-mbc:2+mbc)
    real(kind=8) ::   qr(meqn, 1-mbc:2+mbc)
    real(kind=8) ::   auxl(maux, 1-mbc:2+mbc)
    real(kind=8) ::   auxr(maux, 1-mbc:2+mbc)
    real(kind=8) ::    s_vec(mwaves,         1-mbc:2+mbc)
    real(kind=8) :: wave_vec(meqn,   mwaves, 1-mbc:2+mbc)
    real(kind=8) ::  apdq_vec(meqn,          1-mbc:2+mbc)
    real(kind=8) ::  amdq_vec(meqn,          1-mbc:2+mbc)
    integer :: i, m, mw, mx
    logical, parameter :: debug = .false.

    do i=1,meqn
        qr(i,1) = q_left(i)
        ql(i,2) = q_right(i)
        enddo

    mx = 2

    if (debug) then
        print *, 'qr(:,1):'
        do m=1,meqn
            write(6,60) qr(m,1)
            enddo
        print *, 'ql(i,2):'
        do m=1,meqn
            write(6,60) ql(m,2)
            enddo
        endif

    call rp1(mx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave_vec,s_vec,amdq_vec,apdq_vec)

    ! left and right states in cells 1 and 2, so Riemann solution
    ! is in i=2 element of output vectors:

    do mw=1,mwaves
        s(mw) = s_vec(mw,2)
        do m=1,meqn
            wave(m,mw) = wave_vec(m,mw,2)
            enddo
        enddo
    do m=1,meqn
        amdq(m) = amdq_vec(m,2)
        apdq(m) = apdq_vec(m,2)
        enddo

    if (debug) then
        do i=1-mbc,mx+mbc
            print *, '==> i = ',i
            print *, 's_vec:'
            write(6,60) (s_vec(mw,i), mw=1,mwaves)
60          format(6d16.6)

            print *, 'wave_vec:'
            do m=1,meqn
                write(6,60) (wave_vec(m,mw,i), mw=1,mwaves)
                enddo

            print *, 'amdq:'
            do m=1,meqn
                write(6,60) amdq_vec(m,i)
                enddo
            print *, 'apdq:'
            do m=1,meqn
                write(6,60) apdq_vec(m,i)
                enddo
            enddo
            print *, '==='
        endif

    end subroutine call_rp1
    

end module rp1_driver
