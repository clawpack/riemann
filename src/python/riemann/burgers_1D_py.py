#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for Burgers equation

.. math::
    u_t + \left ( \frac{1}{2} u^2 \right)_x = 0

:Authors:
    Kyle T. Mandli (2009-2-4): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

num_eqn = 1
num_waves = 1

import numpy as np

def burgers_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    Riemann solver for Burgers equation in 1d
         
    *problem_data* should contain -
     - *efix* - (bool) Whether a entropy fix should be used, if not present, 
       false is assumed
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2009-2-4)
    """
        
    num_rp = q_l.shape[1]
    # Output arrays
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.empty( (num_eqn, num_rp) )
    apdq = np.empty( (num_eqn, num_rp) )

    # Basic solve
    wave[0,:,:] = q_r - q_l
    s[0,:] = 0.5 * (q_r[0,:] + q_l[0,:])

    s_index = np.zeros((2,num_rp))
    s_index[0,:] = s[0,:]
    amdq[0,:] = np.min(s_index,axis=0) * wave[0,0,:]
    apdq[0,:] = np.max(s_index,axis=0) * wave[0,0,:]
    
    # Compute entropy fix
    if problem_data['efix']:
        transonic = (q_l[0,:] < 0.0) * (q_r[0,:] > 0.0)
        amdq[0,transonic] = -0.5 * q_l[0,transonic]**2
        apdq[0,transonic] = 0.5 * q_r[0,transonic]**2

    return wave, s, amdq, apdq
