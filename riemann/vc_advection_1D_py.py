#!/usr/bin/env python
# encoding: utf-8
r"""
Variable coefficient 1D advection Riemann solver

.. math::
    q_t + u(x) q_x = 0.

Note that this is the color equation, not the conservative advection equation.
"""

from __future__ import absolute_import
import numpy as np

# Riemann solver constants
num_eqn = 1
num_waves = 1

def vc_advection_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""Basic 1d advection riemann solver
    
    *aux(i)* should contain -
     - *u(x_i)* - (float) advection speed
    
    See :ref:`petclaw_rp` for more details.
    
    :Version: 1.0 (2010-10-10)
    """
    
    # Number of Riemann problems we are solving
    num_rp = q_l.shape[1]
    
    # Return values
    wave = np.empty( (num_eqn,num_waves,num_rp) )
    s = np.empty( (num_waves,num_rp) )
    amdq = np.zeros( (num_eqn,num_rp) )
    apdq = np.zeros( (num_eqn,num_rp) )
    
    wave[0,0,:] = q_r[0,:] - q_l[0,:]

    s[0,:] = aux_l[0,:]

    apdq[0,:] = (aux_l[0,:]>0)*s[0,:] * wave[0,0,:]
    amdq[0,:] = (aux_l[0,:]<0)*s[0,:] * wave[0,0,:]

    return wave, s, amdq, apdq
