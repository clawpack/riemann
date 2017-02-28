
from __future__ import absolute_import
from __future__ import print_function
import test_rp1_meqn2
import numpy

rho = 1.
bulk = 1.
test_rp1_meqn2.cparam.rho = rho
test_rp1_meqn2.cparam.bulk = bulk
test_rp1_meqn2.cparam.cc = numpy.sqrt(bulk/rho)
test_rp1_meqn2.cparam.zz = numpy.sqrt(bulk*rho)

ql = [2.,0.]
qr = [1.,0.]
wave,s,amdq,apdq = test_rp1_meqn2.rp1_driver.call_rp1(2,ql,qr)

print("s = ",s)
print("wave = ")
print(wave)
