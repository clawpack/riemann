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

- Add the following to src/python/riemann/__init__.py:
```
import vc_elasticity_2D
```
- Add 'elasticity' to two_d_riemann in src/python/riemann/setup.py.
- Add appropriate entries for num_eqn and num_waves in src/python/riemann/static.py.
