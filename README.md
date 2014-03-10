# Riemann Solver Repository

This repository is a centralized location for all Clawpack-compatible Riemann
solvers.  If you have developed a Riemann solver that is not already here,
please send it to us or issue a pull request.  The format for Riemann solvers
has changed significantly since Clawpack 4.3, but if you have a 4.3-style solver,
send it along and we will update it.

# Adding a Riemann solver

When adding a new Riemann solver, in addition to adding the normal and
(optionally) transverse solver code files in src/, you should do the following
in order to ensure the new solver is importable in PyClaw.  To make things
concrete, suppose you are committing a 2D elasticity solver, so your Fortran 
files are named rpn2_elasticity.f90 and rpt2_elasticity.f90.  Then you should:

1. Add the following to src/python/riemann/__init__.py:
```
import vc_elasticity_2D
```
2. Add 'vc_elasticity' to two_d_riemann in src/python/riemann/setup.py.
3. Add appropriate entries for num_eqn and num_waves in src/python/riemann/static.py.

If you are adding a Python solver, then you should have a vc_elasticity_2D_py.py 
file in src/python/riemann. Assuming it contains a function vc_elasticity_2D, then 
you should do:

1. Add the following to src/python/riemann/__init__.py:
```
from vc_elasticity_2D_py import vc_elasticity_2D
```
2. Add appropriate entries for num_eqn and num_waves in src/python/riemann/static.py.
