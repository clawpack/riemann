#!/usr/bin/env python
# encoding: utf-8
"""
Wave propagation Riemann solvers implemented in Python and Fortran.
"""

rp_solver_list_1d = []
rp_solver_list_2d = []
rp_solver_list_3d = []

# Import 1d Riemann solvers
from rp_advection import rp_advection_1d
from rp_vc_advection import rp_vc_advection_1d
from rp_acoustics import rp_acoustics_1d
from rp_burgers import rp_burgers_1d
from rp_shallow import rp_shallow_roe_1d, rp_shallow_hll_1d, rp_shallow_exact_1d
from rp_euler import rp_euler_roe_1d
from rp_nonlinear_elasticity import rp_nonlinear_elasticity_1d
import static

try:
    import rp1_acoustics
    import rp1_advection
    import rp1_burgers
    import rp1_euler_with_efix
    import rp1_nonlinear_elasticity_fwave
    import rp1_reactive_euler_with_efix
    import rp1_shallow_roe_with_efix
    import rp1_traffic
    import rp2_acoustics
    import rp2_advection
    import rp2_euler_mapgrid
    import rp2_euler_5wave
    import rp2_euler_4wave
    import rp2_kpp
    import rp2_psystem
    import rp2_shallow_roe_with_efix
    import rp2_shallow_sphere
    import rp2_vc_acoustics
    import rp2_vc_advection
    import rp3_vc_acoustics
except ImportError:
    print 'Warning: Some Riemannn solvers not found.  Remember to make in RIEMANN/src/python/riemann.'
