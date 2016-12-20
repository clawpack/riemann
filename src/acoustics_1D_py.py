#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for constant coefficient acoustics.

.. math:: q_t + A q_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} p(x,t) \\ u(x,t) \end{array} \right ]
  
and the coefficient matrix is 

.. math::

    A = \left [\begin{matrix}
    0      & K\\
    1/\rho & 0 
    \end{matrix} \right ]

The parameters :math:`\rho =` density and :math:`K =` bulk modulus are used
to calculate the impedence :math:`= Z` and speed of sound `= c`.

:Authors:
    Kyle T. Mandli (2009-02-03): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from __future__ import absolute_import
from six.moves import range
num_eqn = 2
num_waves = 2

def acoustics_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    Basic 1d acoustics riemann solver, with interleaved arrays
    
    *problem_data* is expected to contain -
     - *zz* - (float) Impedence
     - *cc* - (float) Speed of sound
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2009-02-03)
    """
    import numpy as np

    # Convenience
    num_rp = np.size(q_l,1)

    # Return values
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.empty( (num_eqn, num_rp) )
    apdq = np.empty( (num_eqn, num_rp) )
    
    # Local values
    delta = np.empty(np.shape(q_l))
    
    delta = q_r - q_l
    a1 = (-delta[0,:] + problem_data['zz']*delta[1,:]) / (2.0 * problem_data['zz'])
    a2 = (delta[0,:] + problem_data['zz']*delta[1,:]) / (2.0 * problem_data['zz'])
        
    # Compute the waves
    # 1-Wave
    wave[0,0,:] = -a1 * problem_data['zz']
    wave[1,0,:] = a1
    s[0,:] = -problem_data['cc']
        
    # 2-Wave
    wave[0,1,:] = a2 * problem_data['zz']
    wave[1,1,:] = a2
    s[1,:] = problem_data['cc']
    
    # Compute the left going and right going fluctuations
    for m in range(num_eqn):
        amdq[m,:] = s[0,:] * wave[m,0,:]
        apdq[m,:] = s[1,:] * wave[m,1,:]
    
    return wave, s, amdq, apdq
