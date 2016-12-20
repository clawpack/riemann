#!/usr/bin/env python
# encoding: utf-8
r"""
F-wave Riemann solver for nonlinear elasticity in heterogeneous media

.. math:: q_t + f(q,x)_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} \epsilon(x,t) \\ \rho(x) u(x,t) \end{array} \right ]
  
and the flux vector is

.. math::

    f(q,x) = \left [ \begin{array}{c} -u \\ \sigma(\epsilon,x) \end{array} \right ]

:Authors:
    David I. Ketcheson (2010-11-06): Initial version
    David I. Ketcheson (2011-05-19): Interleaved
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from __future__ import absolute_import
import numpy as np
from six.moves import range

def nonlinear_elasticity_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    1d nonlinear elasticity riemann solver
    
    *aux* is expected to contain -
     - aux[0,i] - density in cell i
     - aux[1,i] - bulk modulus in cell i
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2010-11-06)
    """
    
    num_eqn = 2
    num_waves = 2
    # Convenience
    num_rp = np.size(q_l,1)

    # Set up arrays for return values
    fwave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.empty( (num_eqn, num_rp) )
    apdq = np.empty( (num_eqn, num_rp) )
    
    #Linearized bulk modulus, sound speed, and impedance:
    bulkl = sigmap(q_l[0,:],aux_l[1,:])
    bulkr = sigmap(q_r[0,:],aux_r[1,:])
    cl = np.sqrt(bulkl/aux_l[0,:])
    cr = np.sqrt(bulkr/aux_r[0,:])
    zl = cl*aux_l[0,:]
    zr = cr*aux_r[0,:]

    #Jumps:
    du   = q_r[1,:]/aux_r[0,:]-q_l[1,:]/aux_l[0,:]
    dsig = sigma(q_r[0,:],aux_r[1,:]) - sigma(q_l[0,:],aux_l[1,:])

    b1 = - (zr*du + dsig) / (zr+zl)
    b2 = - (zl*du - dsig) / (zr+zl)

    # Compute the f-waves
    # 1-Wave
    fwave[0,0,:] = b1
    fwave[1,0,:] = b1 * zl
    s[0,:] = -cl
        
    # 2-Wave
    fwave[0,1,:] = b2
    fwave[1,1,:] = b2 * (-zr)
    s[1,:] = cr
    
    # Compute the left going and right going fluctuations
    for m in range(num_eqn):
        amdq[m,:] = fwave[m,0,:]
        apdq[m,:] = fwave[m,1,:]
    
    return fwave, s, amdq, apdq

def sigma(eps,K):
    return np.exp(K*eps)-1.0

def sigmap(eps,K):
    return K*np.exp(K*eps)
