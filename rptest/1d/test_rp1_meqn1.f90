

program test_rp1

    use rp1_driver, only: call_rp1

    implicit none

    integer, parameter :: meqn = 1
    integer, parameter :: mwaves = 1
    real(kind=8) ::   q_left(meqn)
    real(kind=8) ::   q_right(meqn)
    real(kind=8) ::    s(mwaves)
    real(kind=8) :: wave(meqn,   mwaves)
    real(kind=8) ::  apdq(meqn)
    real(kind=8) ::  amdq(meqn)
    integer :: i, m, mw


    write(6,61) meqn
 61 format('Need to input meqn = ',i3,'  values for each state')
    write(6,*) 'input q_left:'
    read(5,*) (q_left(i), i=1,meqn)
    write(6,*) 'input q_right:'
    read(5,*) (q_right(i), i=1,meqn)

    call call_rp1(meqn,mwaves,q_left,q_right,wave,s,amdq,apdq)



    print *, 'speeds:'
    write(6,60) (s(mw), mw=1,mwaves)
60  format(6d16.6)

    print *, 'waves:'
    do m=1,meqn
        write(6,60) (wave(m,mw), mw=1,mwaves)
        enddo

    print *, 'amdq and apdq:'
    do m=1,meqn
        write(6,60) amdq(m), apdq(m)
        enddo


end program test_rp1
