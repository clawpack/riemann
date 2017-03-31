#!/usr/bin/env python
# encoding: utf-8
"""
Wave propagation Riemann solvers implemented in Python and Fortran.
"""

from __future__ import absolute_import
from __future__ import print_function
rp_solver_list_1d = []
rp_solver_list_2d = []
rp_solver_list_3d = []

# Import 1d Riemann solvers
from . import advection_1D_py
from . import vc_advection_1D_py
from . import acoustics_1D_py
from . import burgers_1D_py
from . import shallow_1D_py
from . import euler_1D_py
from . import nonlinear_elasticity_1D_py
from . import static

from . import acoustics_1D_constants
from . import acoustics_variable_1D_constants
from . import advection_1D_constants
from . import burgers_1D_constants
from . import euler_with_efix_1D_constants
from . import nonlinear_elasticity_fwave_1D_constants
from . import reactive_euler_with_efix_1D_constants
from . import shallow_roe_with_efix_1D_constants
from . import shallow_roe_tracer_1D_constants
from . import traffic_1D_constants
from . import traffic_vc_1D_constants
from . import traffic_vc_tracer_1D_constants
from . import acoustics_2D_constants
from . import acoustics_mapped_2D_constants
from . import advection_2D_constants
from . import burgers_2D_constants
from . import euler_mapgrid_2D_constants
from . import euler_5wave_2D_constants
from . import euler_4wave_2D_constants
from . import kpp_2D_constants
from . import psystem_2D_constants
from . import shallow_roe_with_efix_2D_constants
from . import shallow_sphere_2D_constants
from . import vc_acoustics_2D_constants
from . import vc_advection_2D_constants
from . import vc_elasticity_2D_constants
from . import vc_acoustics_3D_constants
from . import euler_3D_constants
from . import burgers_3D_constants
from . import vc_advection_3D_constants


try:
    from . import acoustics_1D
    from . import acoustics_variable_1D
    from . import acoustics_1D_ptwise
    from . import advection_1D
    from . import advection_1D_ptwise
    from . import burgers_1D
    from . import euler_with_efix_1D
    from . import nonlinear_elasticity_fwave_1D
    from . import reactive_euler_with_efix_1D
    from . import shallow_roe_with_efix_1D
    from . import shallow_bathymetry_fwave_1D
    from . import shallow_roe_tracer_1D
    from . import traffic_1D
    from . import traffic_vc_1D
    from . import traffic_vc_fwave_1D
    from . import traffic_vc_tracer_1D
    from . import acoustics_2D
    from . import acoustics_mapped_2D
    from . import acoustics_2D_ptwise
    from . import advection_2D
    from . import burgers_2D
    from . import euler_mapgrid_2D
    from . import euler_5wave_2D
    from . import euler_4wave_2D
    from . import kpp_2D
    from . import psystem_2D
    from . import shallow_roe_with_efix_2D
    from . import shallow_bathymetry_fwave_2D
    from . import sw_aug_2D
    from . import shallow_sphere_2D
    from . import vc_acoustics_2D
    from . import vc_advection_2D
    from . import vc_elasticity_2D
    from . import vc_acoustics_3D
    from . import euler_3D
    from . import burgers_3D
    from . import vc_advection_3D
except ImportError as e:
    import traceback
    print("********************************************************************")
    print('Warning: Some Riemannn solvers were not able to be imported.')
    print(' Did you run "pip install" in your clawpack directory?')
    traceback.print_exc()
    print("********************************************************************")

import os
if os.path.exists('./layered_shallow_water_1D.so'):
    import layered_shallow_water_1D
    from . import layered_shallow_water_1D_constants
