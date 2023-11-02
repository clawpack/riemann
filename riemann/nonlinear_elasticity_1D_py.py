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

Two stress-strain relations are included.  The default is

..  math:: \sigma(\epsilon,K) = \exp(K(x) \epsilon) - 1

By setting problem_data['stress_relation'] = 'quadratic' one can use instead

..  math:: \sigma(\epsilon,K_1,K_2) = K_1 \epsilon + K_2 \epsilon^2

For the latter relation, K_2 is stored in a third aux entry (aux[2,:]).
"""

import numpy as np

strain = 0
momentum = 1
density = 0
K = 1
K1 = 1
K2 = 2


def nonlinear_elasticity_1D(q_l, q_r, aux_l, aux_r, problem_data):
    r"""
    1d nonlinear elasticity riemann solver

    *aux* is expected to contain -
     - aux[0,i] - density in cell i
     - aux[1,i] - bulk modulus in cell i
     - aux[2,i] - K2 (needed for quadratic stress relation only)

    See :ref:`pyclaw_rp` for more details.
    """

    num_eqn = 2
    num_waves = 2
    # Convenience
    num_rp = np.size(q_l, 1)

    # Set up arrays for return values
    fwave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.empty( (num_eqn, num_rp) )
    apdq = np.empty( (num_eqn, num_rp) )

    sigma_l = sigma(q_l, aux_l, problem_data)
    sigma_r = sigma(q_r, aux_r, problem_data)
    bulk_l = sigmap(q_l, aux_l, problem_data)
    bulk_r = sigmap(q_r, aux_r, problem_data)

    # Linearized bulk modulus, sound speed, and impedance:
    cl = np.sqrt(bulk_l/aux_l[density,:])
    cr = np.sqrt(bulk_r/aux_r[density,:])
    zl = cl*aux_l[density,:]
    zr = cr*aux_r[density,:]

    # Jumps:
    du = q_r[momentum,:]/aux_r[density,:]-q_l[momentum,:]/aux_l[density,:]
    dsig = sigma_r - sigma_l

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


def sigma(q, aux, problem_data):
    r"""Stress-strain relation."""
    stress_relation = problem_data.get('stress_relation', 'exponential')
    if stress_relation == 'exponential':
        return np.exp(aux[K,:]*q[strain,:]) - 1.
    elif stress_relation == 'quadratic':
        return aux[K1,:]*q[strain,:] + aux[K2,:]*q[strain,:]**2


def sigmap(q, aux, problem_data):
    r"""Derivative of stress-strain relation w.r.t. strain."""
    stress_relation = problem_data.get('stress_relation', 'exponential')
    if stress_relation == 'exponential':
        return aux[K,:]*np.exp(aux[K,:]*q[strain,:])
    elif stress_relation == 'quadratic':
        return aux[K1,:] + 2*aux[K2,:]*q[strain,:]
