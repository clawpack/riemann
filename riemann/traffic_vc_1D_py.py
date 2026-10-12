#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solver for variable coefficient LWR traffic model in one dimension.

.. math::
    q_t + f(q,x)_x = 0

with the flux function

.. math::

    f(q,x) = u(x) q (1-q)
"""

import numpy as np

def traffic_vc_1D_py(q_l, q_r, aux_l, aux_r, problem_data) :
    """
    Godunov scheme for solving the 1D variable-coefficient traffic problem.

    Variables `aux_l` and `aux_r` contain speeds in neighbouring cells.

    :Version: 1.0 (28/01/2024) by Kaloyan Danovski
    """
    # Number of Riemann problems we are solving
    num_rp = q_l.shape[1]

    # Return values
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Extract velocities
    v_l, v_r = aux_l[0,:], aux_r[0,:]
    
    # Calculate fluxes
    f_l = v_l * q_l * (1 - q_l)
    f_r = v_r * q_r * (1 - q_r)
    
    # Calculate characteristic speeds
    c_l = v_l * (1 - 2 * q_l)
    c_r = v_r * (1 - 2 * q_r)

    # Calculate wave and speeds
    wave[0,0,:] = q_r[0,:] - q_l[0,:]
    s[0,:] = 0.5 * (c_l + c_r)

    # Define condition masks
    # (vector-programming logic)
    cr_pos = c_r > 0
    cr_neg = ~cr_pos
    cl_pos = c_l > 0
    cl_neg = ~cl_pos
    f_lggr = f_l > f_r
    f_lllr = ~f_lggr
    f_lggvr = f_l > v_r/4
    f_rggvl = f_r > v_l/4
    v_lggr = v_l >= v_r
    v_lllr = ~v_lggr

    transonic = cl_neg * cr_pos

    # Compute middle state
    f0 = np.zeros_like(f_r)
    f0 += cl_pos * cr_pos *  f_lggvr          * v_r/4       # Left-going shock, right-going raref
    f0 += cl_neg * cr_neg *  f_rggvl          * v_l/4       # Right-going shock, left-going raref
    f0 +=          cr_neg * ~f_lggvr * f_lggr * f_r         # Left-going shock
    f0 += cl_pos *          ~f_rggvl * f_lllr * f_l         # Right-going shock
    f0 += cl_neg * cr_neg * ~f_rggvl * f_lllr * f_r         # Left-going raref
    f0 += cl_pos * cr_pos * ~f_lggvr * f_lggr * f_l         # Right-going raref
    f0 += transonic * v_lggr                  * v_r/4       # Transonic (left)
    f0 += transonic * v_lllr                  * v_l/4       # Transonic (right)

    # Update fluctuations
    amdq[0,:] = f0 - f_l
    apdq[0,:] = f_r - f0

    return wave, s, amdq, apdq