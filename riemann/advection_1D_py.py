#!/usr/bin/env python
# encoding: utf-8
r"""
Simple advection Riemann solvers

Basic advection Riemann solvers of the form (1d)

.. math::
    q_t + A q_x = 0.
    
:Authors:
    Kyle T. Mandli (2008-2-20): Initial version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

# Riemann solver constants
num_eqn = 1
num_waves = 1

import numpy as np

def advection_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""Basic 1d advection riemann solver
    
    *problem_data* should contain -
     - *u* - (float) Determines advection speed
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2008-2-20)
    """
    
    # Number of Riemann problems we are solving
    num_rp = q_l.shape[1]

    # Return values
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )
 
    wave[0,0,:] = q_r[0,:] - q_l[0,:]
    s[0,:] = problem_data['u']
    if problem_data['u'] > 0:
        apdq[0,:] = s[0,:] * wave[0,0,:]
    else:
        amdq[0,:] = s[0,:] * wave[0,0,:]

    return wave, s, amdq, apdq
