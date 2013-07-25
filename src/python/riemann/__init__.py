#!/usr/bin/env python
# encoding: utf-8
"""
Wave propagation Riemann solvers implemented in Python and Fortran.
"""

rp_solver_list_1d = []
rp_solver_list_2d = []
rp_solver_list_3d = []

# Import 1d Riemann solvers
from advection_1D_py import advection_1D
from vc_advection_1D_py import vc_advection_1D
from acoustics_1D_py import acoustics_1D
from burgers_1D_py import burgers_1D
from shallow_1D_py import shallow_roe_1D, shallow_hll_1D, shallow_exact_1D
from euler_1D_py import euler_roe_1D
from nonlinear_elasticity_1D_py import nonlinear_elasticity_1D
import static

try:
    import acoustics_1D
    import advection_1D
    import burgers_1D
    import euler_with_efix_1D
    import nonlinear_elasticity_fwave_1D
    import reactive_euler_with_efix_1D
    import shallow_roe_with_efix_1D
    import layered_shallow_water_1D
    import traffic_1D
    import acoustics_2D
    import advection_2D
    import burgers_2D
    import euler_mapgrid_2D
    import euler_5wave_2D
    import euler_4wave_2D
    import kpp_2D
    import psystem_2D
    import shallow_roe_with_efix_2D
    import shallow_sphere_2D
    import vc_acoustics_2D
    import vc_advection_2D
    import vc_acoustics_3D
    import euler_3D
    import burgers_3D
    import vc_advection_3D
except ImportError as e:
    import traceback
    print "********************************************************************"
    print 'Warning: Some Riemannn solvers were not able to be imported.'
    print '         Did you run "pip install" in your clawpack directory?'
    traceback.print_exc()
    print "********************************************************************"
