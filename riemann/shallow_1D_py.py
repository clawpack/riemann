#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for the shallow water equations.

The available solvers are:
 * Roe - Use Roe averages to caluclate the solution to the Riemann problem
 * HLL - Use a HLL solver
 * Exact - Use a newton iteration to calculate the exact solution to the
        Riemann problem

.. math::
    q_t + f(q)_x = 0

where

.. math::
    q(x,t) = \left [ \begin{array}{c} h \\ h u \end{array} \right ],

the flux function is

.. math::
    f(q) = \left [ \begin{array}{c} h u \\ hu^2 + 1/2 g h^2 \end{array}\right ].

and :math:`h` is the water column height, :math:`u` the velocity and :math:`g`
is the gravitational acceleration.
"""

import numpy as np

num_eqn = 2
num_waves = 2


def shallow_roe_1D(q_l, q_r, aux_l, aux_r, problem_data):
    r"""
    Roe shallow water solver in 1d::

        ubar = (sqrt(u_l) + sqrt(u_r)) / (sqrt(h_l) + sqrt(h_r))
        cbar = sqrt( 0.5 * g * (h_l + h_r))

        W_1 = |      1      |  s_1 = ubar - cbar
              | ubar - cbar |

        W_2 = |      1      |  s_1 = ubar + cbar
              | ubar + cbar |

        a1 = 0.5 * ( - delta_hu + (ubar + cbar) * delta_h ) / cbar
        a2 = 0.5 * (   delta_hu - (ubar - cbar) * delta_h ) / cbar

    *problem_data* should contain:
     - *g* - (float) Gravitational constant
     - *efix* - (bool) Boolean as to whether a entropy fix should be used, if
       not present, false is assumed

    :Version: 1.0 (2009-02-05)
    """

    # Array shapes
    num_rp = q_l.shape[1]

    # Output arrays
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.zeros( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Compute roe-averaged quantities
    ubar = ( (q_l[1,:]/np.sqrt(q_l[0,:]) + q_r[1,:]/np.sqrt(q_r[0,:])) /
             (np.sqrt(q_l[0,:]) + np.sqrt(q_r[0,:])) )
    cbar = np.sqrt(0.5 * problem_data['grav'] * (q_l[0,:] + q_r[0,:]))

    # Compute Flux structure
    delta = q_r - q_l
    a1 = 0.5 * (-delta[1,:] + (ubar + cbar) * delta[0,:]) / cbar
    a2 = 0.5 * ( delta[1,:] - (ubar - cbar) * delta[0,:]) / cbar

    # Compute each family of waves
    wave[0,0,:] = a1
    wave[1,0,:] = a1 * (ubar - cbar)
    s[0,:] = ubar - cbar

    wave[0,1,:] = a2
    wave[1,1,:] = a2 * (ubar + cbar)
    s[1,:] = ubar + cbar

    if problem_data['efix']:
        raise NotImplementedError("Entropy fix has not been implemented.")
    else:
        s_index = np.zeros((2,num_rp))
        for m in range(num_eqn):
            for mw in range(num_waves):
                s_index[0,:] = s[mw,:]
                amdq[m,:] += np.min(s_index,axis=0) * wave[m,mw,:]
                apdq[m,:] += np.max(s_index,axis=0) * wave[m,mw,:]

    return wave, s, amdq, apdq

def shallow_hll_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    HLL shallow water solver ::


        W_1 = Q_hat - Q_l    s_1 = min(u_l-c_l,u_l+c_l,lambda_roe_1,lambda_roe_2)
        W_2 = Q_r - Q_hat    s_2 = max(u_r-c_r,u_r+c_r,lambda_roe_1,lambda_roe_2)

        Q_hat = ( f(q_r) - f(q_l) - s_2 * q_r + s_1 * q_l ) / (s_1 - s_2)

    *problem_data* should contain:
     - *g* - (float) Gravitational constant

    :Version: 1.0 (2009-02-05)
    """
    # Array shapes
    num_rp = q_l.shape[1]
    num_eqn = 2
    num_waves = 2

    # Output arrays
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Compute Roe and right and left speeds
    ubar = ( (q_l[1,:]/np.sqrt(q_l[0,:]) + q_r[1,:]/np.sqrt(q_r[0,:])) /
        (np.sqrt(q_l[0,:]) + np.sqrt(q_r[0,:])) )
    cbar = np.sqrt(0.5 * problem_data['grav'] * (q_l[0,:] + q_r[0,:]))
    u_r = q_r[1,:] / q_r[0,:]
    c_r = np.sqrt(problem_data['grav'] * q_r[0,:])
    u_l = q_l[1,:] / q_l[0,:]
    c_l = np.sqrt(problem_data['grav'] * q_l[0,:])

    # Compute Einfeldt speeds
    s_index = np.empty((4,num_rp))
    s_index[0,:] = ubar+cbar
    s_index[1,:] = ubar-cbar
    s_index[2,:] = u_l + c_l
    s_index[3,:] = u_l - c_l
    s[0,:] = np.min(s_index,axis=0)
    s_index[2,:] = u_r + c_r
    s_index[3,:] = u_r - c_r
    s[1,:] = np.max(s_index,axis=0)

    # Compute middle state
    q_hat = np.empty((2,num_rp))
    q_hat[0,:] = ((q_r[1,:] - q_l[1,:] - s[1,:] * q_r[0,:]
                            + s[0,:] * q_l[0,:]) / (s[0,:] - s[1,:]))
    q_hat[1,:] = ((q_r[1,:]**2/q_r[0,:] + 0.5 * problem_data['grav'] * q_r[0,:]**2
                - (q_l[1,:]**2/q_l[0,:] + 0.5 * problem_data['grav'] * q_l[0,:]**2)
                - s[1,:] * q_r[1,:] + s[0,:] * q_l[1,:]) / (s[0,:] - s[1,:]))

    # Compute each family of waves
    wave[:,0,:] = q_hat - q_l
    wave[:,1,:] = q_r - q_hat

    # Compute variations
    s_index = np.zeros((2,num_rp))
    for m in range(num_eqn):
        for mw in range(num_waves):
            s_index[0,:] = s[mw,:]
            amdq[m,:] += np.min(s_index,axis=0) * wave[m,mw,:]
            apdq[m,:] += np.max(s_index,axis=0) * wave[m,mw,:]

    return wave, s, amdq, apdq


def shallow_fwave_1d(q_l, q_r, aux_l, aux_r, problem_data):
    r"""Shallow water Riemann solver using fwaves

    Also includes support for bathymetry but be wary if you think you might have
    dry states as this has not been tested.

    *problem_data* should contain:
     - *grav* - (float) Gravitational constant
     - *dry_tolerance* - (float) Set velocities to zero if h is below this
       tolerance.
     - *sea_level* - (float) Datum from which the dry-state is calculated.

    :Version: 1.0 (2014-09-05)
    :Version: 2.0 (2017-03-07)
    """

    g = problem_data['grav']
    dry_tolerance = problem_data['dry_tolerance']
    sea_level = problem_data['sea_level']

    num_rp = q_l.shape[1]
    num_eqn = 2
    num_waves = 2

    # Output arrays
    fwave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Extract state
    u_l = np.where(q_l[0, :] > dry_tolerance,
                   q_l[1, :] / q_l[0, :], 0.0)
    u_r = np.where(q_r[0, :] > dry_tolerance,
                   q_r[1, :] / q_r[0, :], 0.0)
    phi_l = q_l[0, :] * u_l**2 + 0.5 * g * q_l[0, :]**2
    phi_r = q_r[0, :] * u_r**2 + 0.5 * g * q_r[0, :]**2
    h_bar = 0.5 * (q_l[0, :] + q_r[0, :])

    # Speeds
    u_hat = (np.sqrt(g * q_l[0, :]) * u_l + np.sqrt(g * q_r[0, :]) * u_r)      \
            / (np.sqrt(g * q_l[0, :]) + np.sqrt(g * q_r[0, :]))
    c_hat = np.sqrt(g * h_bar)
    s[0, :] = np.amin(np.vstack((u_l - np.sqrt(g * q_l[0, :]),
                                 u_hat - c_hat)), axis=0)
    s[1, :] = np.amax(np.vstack((u_r + np.sqrt(g * q_r[0, :]),
                                 u_hat + c_hat)), axis=0)

    delta1 = q_r[1, :] - q_l[1, :]
    delta2 = phi_r - phi_l + g * h_bar * (aux_r[0, :] - aux_l[0, :])

    beta1 = (s[1, :] * delta1 - delta2) / (s[1, :] - s[0, :])
    beta2 = (delta2 - s[0, :] * delta1) / (s[1, :] - s[0, :])

    fwave[0, 0, :] = beta1
    fwave[1, 0, :] = beta1 * s[0, :]
    fwave[0, 1, :] = beta2
    fwave[1, 1, :] = beta2 * s[1, :]

    for m in range(num_eqn):
        for mw in range(num_waves):
            amdq[m, :] += (s[mw, :] < 0.0) * fwave[m, mw, :]
            apdq[m, :] += (s[mw, :] > 0.0) * fwave[m, mw, :]

            amdq[m, :] += (s[mw, :] == 0.0) * fwave[m, mw, :] * 0.5
            apdq[m, :] += (s[mw, :] == 0.0) * fwave[m, mw, :] * 0.5

    return fwave, s, amdq, apdq


def shallow_exact_1D(q_l, q_r, aux_l, aux_r, problem_data):
    r"""
    Exact shallow water Riemann solver

    .. warning::
        This solver has not been implemented.

    """
    raise NotImplementedError("The exact swe solver has not been implemented.")
