#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solver for Hughes' model of pedestrian motion in one dimension.

The model describes crowd dynamics in an evacuation scenario, i.e. wish a shared goal of escaping a given space. In its simplest form, the model reads

.. math::
    q_t - (\text{sgn}(\phi_x) f(q))_x = 0

where
* :math:`f(q)` is the model flux defined as :math:`f(q) = qv(q) = q(1-q)`,
* :math:`phi` is a potential function defined by the following eikonal equation.

.. math::
    |\phi_x| = c(q) = 1/v(q)

This is determined with the so-called Fast Sweeping Algorithm, which however requires sequentially sweeping the domain and is thus slow in Python. A better, more "naive" alternative is provided in this implementation.

The point where the derivative of the potential :math:`phi_x` changes sign is called the *turning point* $\xi(t)$. A more general formulation of the model reads

.. math::
    q_t + F(t,x,q)_x = 0

where :math:`F(t,x,q) = (x-\xi(t))f(q)`. This problem is a discontinuous, space-time dependent flux with a velocity given by a non-stationary turning point. The velocity is equal to :math:`\text{sgn}(\phi_x)=\pm1`.

To solve it, at every time step we have to:
1. Given :math:`q(t,x)`, solve the eikonal equation to obtain :math:`\phi`.
2. Given :math:`\phi` (i.e. knowing :math:`\xi` and thus the velocity), evolve the solution by solving the Riemann problem with a suitable scheme.

The Godunov and Local Lax-Friedrichs methods are implemented.
"""

import numpy as np

# Define useful functions
v = lambda q : 1 - q        # Velocity
c = lambda q : 1 / v(q)     # Cost function
f = lambda q : q * v(q)     # Flux
df = lambda q : 1 - 2 * q   # Characteristic speed

def sweep(q, phi, js, dx, err_th=10e-10) :
    r"""
    Fast Sweeping Algorithm for solving an eikonal differential equation in 1D.
    
    In this case the equation reads:
    
    .. math::
        |\phi_x| = c(q) = 1/v(q)

    The algorithm proceeds by sequentially sweeping the domain in both directions, defined by the :code:`js` parameter. (Logic for updating the potential is available below.)

    It stops when a fixed error threshold :code:`err_th` is met for :math:`|\phi^\text{old} - \phi^\text{new}|`.

    .. note::
        Updates phi in-place.
    """
    err = 0
    n_sweeps = 1    # Used for performance assessment; means that program returned after the n-th sweep

    while True :    # Will break out when error is met
        for j in js :
            phi_old = phi[j]
            phi_xmin = min(phi[j-1], phi[j+1])
            phi_sol = (1 / (1-q[j])) * dx + phi_xmin
            phi[j] = min(phi_old, phi_sol)   # Update phi
            
            err += phi_old - phi[j]
        if err < err_th : return n_sweeps   # Exit when threshold is met
        # Otherwise...
        n_sweeps += 1       # increment,
        err = 0             # reset error,
        js = js[::-1]       # and set up to sweep from the other side.

def riemann_solve_pos(q_l, q_r, llf=False) :
    '''
    Riemann solver for the half-problem (after splitting around turning point) with constant speed = 1.
    
    Code assumes the left half-problem is transposed and reversed such that the first index in the arguments stands for the boundary of the cell containing the turning point :math:`\xi`. That is, :code:`q_l[0]` is the flux in the turning cell.

    To properly account for the flux going in the opposite direction from the other boundary of the turning point :math:`\xi`, we need to add the corresponding fluctuations as a left-going fluctuation in the first boundary, i.e. :code:`amdq[0]`.

    The method used (Godunov or Local Lax-Friedrichs) is given by the boolean parameter :code:`llf`, which is :code:`False` by default.
    '''
    # Calculate fluxes
    f_l, f_r = f(q_l), f(q_r)
    # Calculate characteristic speeds
    c_l, c_r = df(q_l), df(q_r)

    # Define logic masks
    cr_pos = c_r > 0
    cr_neg = ~cr_pos
    cl_pos = c_l > 0
    cl_neg = ~cl_pos
    f_lggr = f_l > f_r
    f_lllr = ~f_lggr
    transonic = cl_neg * cr_pos

    # Compute middle state
    f0 = np.zeros_like(f_r)
    f0 +=          cr_neg * f_lggr * f_r         # Left-going shock
    f0 += cl_pos *          f_lllr * f_l         # Right-going shock
    f0 += cl_neg * cr_neg * f_lllr * f_r         # Left-going raref
    f0 += cl_pos * cr_pos * f_lggr * f_l         # Right-going raref
    f0 += transonic                * 1/4         # Transonic raref

    # Compute fluctuations
    if llf :    # Use Local Lax-Friedrichs / Rusanov method
        a = np.max(np.stack([np.abs(c_l), np.abs(c_r)]), axis=0)
        q_diff = q_r - q_l
        amdq = f0 - f_l - 0.5 * a * q_diff
        apdq = f_r - f0 + 0.5 * a * q_diff
    else :  # Use Godunov method
        amdq = f0 - f_l
        apdq = f_r - f0

    # Fix for fluctuations, to take into account flux difference
    # inside the cell containing the turning point.
    amdq[0] += f_l[0]

    return amdq, apdq

def hughes_1D_py(q_l, q_r, aux_l, aux_r, problem_data) :
    r"""
    Riemann solver (including solution of eikonal equation for potential) for Hughes' pedestrian model, with a discontinuous, space-time dependent flux.

    The expected :code:`problem_data` parameters are:
    * 'llf' : (bool) Whether to use the Local Lax-Friedrichs method. If :code:`False`, use the classic Godunov method instead.
    * 'sweep' : (bool) Whether to use the Fast Sweeping Algorithm to solve the eikonal equation. If :code:`False`, use the faster, "naive" approach instead.
    * 'N' : (int) Number of cells for spatial discretization (used to calculate dx).

    The problem is solved by splitting the domain (in this case, boundary arrays) into two around the turning point, the left and right half problems, and then solving each separately via a simple traffic model-like Riemann solver with a small fix around the turning point :math:`\xi`.
    """

    # Problem parameters
    num_rp = q_l.shape[1]

    # Return values
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )    # A minus delta q
    apdq = np.zeros( (num_eqn, num_rp) )    # A plus delta q

    ### Compute turning point
    # The variable `k` contains the sign of phi and is used
    # afterwards for splitting the problem
    if problem_data['sweep'] :  # Via fast sweeping algorithm
        q = np.concatenate([q_l[0,:], np.zeros(1)])
        num_ghost = 2
        j_interior = list(range(num_ghost, q.size-num_ghost))
        dx = 2/problem_data['N']

        # Initialize phi
        phi = np.zeros_like(q)
        phi[j_interior] = 100

        # Run algorithm
        n_sweeps = sweep(q, phi, j_interior, dx)

        # Calculate derivative
        sgn_phi = np.sign(phi[1:] - phi[:-1])

        # Compute `k` — sign of derivative
        k = sgn_phi
        k[0], k[-1] = 1, -1    # Small fix, to make splitting logic easier

    else :  # Via naive algorithm
        c_all = 1 / (1 - (q_l[0,:] - 10e-15))  # Will break down if q=1, so we subtract a very small number
        c_total = np.sum(c_all)
        j = 1       # Start sum from 1 to avoid problems in case the total array size is zero
        csum = c_all[0]
        while csum + c_all[j] < c_total/2 :
            csum += c_all[j]
            j += 1  # Index of right boundary of cell with turning point, i.e. xi is located in cell j
        
        xi = j - 1  # Switch to index of left boundary
    
        # Compute `k` — sign of derivative
        k_idx = np.arange(q_l.size)
        k = (k_idx <= xi) * 1 + (k_idx > xi) * -1
    
    # Define masks for splitting problem into left and right of xi
    l_mask = k > 0
    r_mask = k < 0

    # Compute waves
    wave[0,0,:] = q_r[0,:] - q_l[0,:]
    # Compute characteristic speeds
    q_factor = 1 - q_l[0,:] + q_r[0,:]
    s[0,:] = r_mask * q_factor + l_mask * -q_factor
    
    # Solve via problem splitting
    # At the left (where real velocity = -1), we transpose boundary directions and reverse arrays
    amdq_left, apdq_left = riemann_solve_pos(q_r[0,l_mask][::-1], q_l[0,l_mask][::-1], problem_data['llf'])
    amdq_right, apdq_right = riemann_solve_pos(q_l[0,r_mask], q_r[0,r_mask], problem_data['llf'])
    amdq[0,:] = np.concatenate([apdq_left[::-1], amdq_right])   # Need to reverse again when recombining
    apdq[0,:] = np.concatenate([amdq_left[::-1], apdq_right])

    return wave, s, amdq, apdq