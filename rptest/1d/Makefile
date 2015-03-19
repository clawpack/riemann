
#FFLAGS = -W -Wall -fbounds-check -pedantic-errors -std=f2003 -Wunderflow -O -ffpe-trap=zero,overflow,underflow,invalid  -g

FFLAGS =

test1:
	gfortran  rp1_driver.f90 test_rp1_meqn1.f90 ../../src/rp1_burgers.f90 
	./a.out
    
test2:
	gfortran $(FFLAGS) rp1_driver.f90 test_rp1_meqn2.f90 ../../src/rp1_acoustics.f90 
	./a.out
    
test1_f2py:
	f2py -c rp1_driver.f90  ../../src/rp1_burgers.f90 -m test_rp1_meqn1
	@echo "You should now be able to..."
	@echo "    >>> from test_rp1_meqn1 import rp1_driver"
	@echo "    >>> ql = 3; qr = 1"
	@echo "    >>> wave,s,amdq,apdq = rp1_driver.call_rp1(1,ql,qr)"

    
test2_f2py:
	f2py -c rp1_driver.f90  ../../src/rp1_acoustics.f90 -m test_rp1_meqn2
	@echo "You should now be able to..."
	@echo "    python test_rp1_meqn2.py"

clean:
	rm *.o *.mod a.out
