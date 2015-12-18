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

:Authors:
    Kyle T. Mandli (2009-02-05): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

num_eqn = 2
num_waves = 2

def shallow_roe_1D(q_l,q_r,aux_l,aux_r,problem_data):
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
            
    s_index = np.zeros((2,num_rp))
    for m in xrange(num_eqn):
        for mw in xrange(num_waves):
            s_index[0,:] = s[mw,:]
            amdq[m,:] += np.min(s_index,axis=0) * wave[m,mw,:]
            apdq[m,:] += np.max(s_index,axis=0) * wave[m,mw,:]

    if problem_data['efix']:
        # Compute eigenvalues of f'(q)
        def lamb(i,q):
            if (i == 1):
                return q[1,:]/q[0,:] - np.sqrt(problem_data['grav']*q[0,:])
            else:
                return q[1,:]/q[0,:] + np.sqrt(problem_data['grav']*q[0,:])

        # Compute intermediate state
        q_m = q_l + wave[:,0,:]

        # Check for transonic rarefactions
        transonic_1   = (lamb(1,q_l) < 0.0) * (lamb(1,q_m) > 0.0)
        transonic_2   = (lamb(2,q_m) < 0.0) * (lamb(2,q_r) > 0.0)

        # Harten-Hyman entropy fix parameter
        beta = np.zeros(num_rp)
        beta[transonic_1] = (lamb(1,q_m[:,transonic_1]) - s[0,transonic_1]) \
            / (lamb(1,q_m[:,transonic_1]) - lamb(1,q_l[:,transonic_1]))
        beta[transonic_2] = (lamb(2,q_r[:,transonic_2]) - s[1,transonic_2]) \
            / (lamb(2,q_r[:,transonic_2]) - lamb(2,q_m[:,transonic_2]))

        # Update fluctuations
        amdq[:,transonic_1] += beta[transonic_1] * lamb(1,q_l[:,transonic_1]) \
            * wave[:,0,transonic_1] * (s[0,transonic_1] >= 0.0) \
            + (beta[transonic_1] * lamb(1,q_l[:,transonic_1]) \
            - s[0,transonic_1]) * wave[:,0,transonic_1] \
            * (s[0,transonic_1] < 0.0)
        amdq[:,transonic_2] += beta[transonic_2] * lamb(2,q_m[:,transonic_2]) \
            * wave[:,1,transonic_2] * (s[1,transonic_2] >= 0.0) \
            + (beta[transonic_2] * lamb(2,q_m[:,transonic_2]) \
            - s[1,transonic_2]) * wave[:,1,transonic_2] \
            * (s[1,transonic_2] < 0.0)
        apdq[:,transonic_1] += (1. - beta[transonic_1]) \
            * lamb(1,q_m[:,transonic_1]) * wave[:,0,transonic_1] \
            * (s[0,transonic_1] < 0.0) + ((1. - beta[transonic_1]) \
            * lamb(1,q_m[:,transonic_1]) - s[0,transonic_1]) \
            * wave[:,0,transonic_1] * (s[0,transonic_1] >= 0.0)
        apdq[:,transonic_2] += (1. - beta[transonic_2]) \
            * lamb(2,q_r[:,transonic_2]) * wave[:,1,transonic_2] \
            * (s[1,transonic_2] < 0.0) + ((1. - beta[transonic_2]) \
            * lamb(2,q_r[:,transonic_2]) - s[1,transonic_2]) \
            * wave[:,1,transonic_2] * (s[1,transonic_2] >= 0.0)
            
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
    q_hat[0,:] = ((q_r[1,:] - q_l[1,:] - s[1,:] * q_r[0,:] \
                   + s[0,:] * q_l[0,:]) / (s[0,:] - s[1,:]))
    q_hat[1,:] = ((q_r[1,:]**2/q_r[0,:] + 0.5 * problem_data['grav'] * q_r[0,:]**2 \
                   - (q_l[1,:]**2/q_l[0,:] + 0.5 * problem_data['grav'] * q_l[0,:]**2) \
        - s[1,:] * q_r[1,:] + s[0,:] * q_l[1,:]) / (s[0,:] - s[1,:]))

    # Compute each family of waves
    wave[:,0,:] = q_hat - q_l
    wave[:,1,:] = q_r - q_hat
    
    # Compute variations
    s_index = np.zeros((2,num_rp))
    for m in xrange(num_eqn):
        for mw in xrange(num_waves):
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
     - *sea_level* - (float) Datum from which the dry-state is calculated.
            
    :Version: 1.0 (2014-09-05)
    """

    g = problem_data['grav']

    num_rp = q_l.shape[1]
    num_eqn = 2
    num_waves = 2
    
    # Output arrays
    fwave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Extract state
    u_l = np.where(q_l[0,:] - problem_data['sea_level'] > 1e-3, 
                   q_l[1,:] / q_l[0,:], 0.0)
    u_r = np.where(q_r[0,:] - problem_data['sea_level'] > 1e-3, 
                   q_r[1,:] / q_r[0,:], 0.0)
    phi_l = q_l[0,:] * u_l**2 + 0.5 * g * q_l[0,:]**2
    phi_r = q_r[0,:] * u_r**2 + 0.5 * g * q_r[0,:]**2

    # Speeds
    s[0,:] = u_l - np.sqrt(g * q_l[0,:])
    s[1,:] = u_r + np.sqrt(g * q_r[0,:])

    delta1 = q_r[1,:] - q_l[1,:]
    delta2 = phi_r - phi_l + g * 0.5 * (q_r[0,:] + q_l[0,:]) * (aux_r[0,:] - aux_l[0,:])

    beta1 = (s[1,:] * delta1 - delta2) / (s[1,:] - s[0,:])
    beta2 = (delta2 - s[0,:] * delta1) / (s[1,:] - s[0,:])

    fwave[0,0,:] = beta1
    fwave[1,0,:] = beta1 * s[0,:]
    fwave[0,1,:] = beta2
    fwave[1,1,:] = beta2 * s[1,:]

    for m in xrange(num_eqn):
        for mw in xrange(num_waves):
            amdq[m,:] += (s[mw,:] < 0.0) * fwave[m,mw,:]
            apdq[m,:] += (s[mw,:] >= 0.0) * fwave[m,mw,:]

    return fwave, s, amdq, apdq

    
def shallow_exact_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    Exact shallow water Riemann solver

    """
    num_eq = 2
    num_waves = 2
    
    # Parameters
    g = problem_data['grav']

    # Array shapes
    num_rp = q_l.shape[1]

    # Output arrays
    wave = np.zeros( (num_eqn, num_waves, num_rp) )
    s = np.zeros( (num_waves, num_rp) )
    sm = np.zeros( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Set heights and velocities
    h_l, h_r = q_l[0,:], q_r[0,:]
    u_l, u_r = q_l[1,:] / q_l[0,:], q_r[1,:] / q_r[0,:]

    # Set intermediate states
    h_m, u_m = np.zeros(num_rp), np.zeros(num_rp) 

    # Functions defined in George 2008 (Appendix B)
    def phi(x, h_p):
        if (x <= h_p):
            return 2.*(np.sqrt(g*x) - np.sqrt(g*h_p))
        else:
            return (x - h_p)*np.sqrt(0.5*g*(1./x + 1./h_p))

    def psi(x, h_l, h_r, u_l, u_r):
        return phi(x, h_r) + phi(x, h_l) + u_r - u_l

    psi_min, psi_max = np.zeros(num_rp), np.zeros(num_rp)

    # Newton solve to find intermediate state q_m
    for i in xrange(num_rp):
        h_m[i] = scipy.optimize.newton(psi, 1.e-3, \
                                       args=(h_l[i],h_r[i],u_l[i],u_r[i]))
        u_m[i] = (u_l[i] - phi(h_m[i], h_l[i]))
        h_min, h_max = min(h_l[i], h_r[i]), max(h_l[i], h_r[i])
        psi_min[i] = psi(h_min, h_l[i], h_r[i], u_l[i], u_r[i])
        psi_max[i] = psi(h_max, h_l[i], h_r[i], u_l[i], u_r[i])

    # Compute Roe and right and left speeds
    ubar = ( (q_l[1,:]/np.sqrt(q_l[0,:]) + q_r[1,:]/np.sqrt(q_r[0,:]))
             / (np.sqrt(q_l[0,:]) + np.sqrt(q_r[0,:])) )
    cbar = np.sqrt(0.5*g*(q_l[0,:] + q_r[0,:]))
    u_r = q_r[1,:]/q_r[0,:]
    c_r = np.sqrt(g*q_r[0,:])
    u_l = q_l[1,:]/q_l[0,:]
    c_l = np.sqrt(g*q_l[0,:])

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
    
    # Determine characteristic structure for each Riemann problem       
    all_shock = (psi_min <= psi_max)*(psi_max <= 0.0)
    one_rar = (psi_min < 0.0)*(psi_max >= 0.0)*(h_l > h_r)
    two_rar = (psi_min < 0.0)*(psi_max > 0.0)*(h_l < h_r)
    all_rar = (0.0 <= psi_min)*(psi_min < psi_max)    

    # qt1 and qt2 are transonic rarefactions in the 1- and 2-wave, respectively.   
    qt1, qt2 = np.zeros( (num_eqn, num_rp) ), np.zeros( (num_eqn, num_rp) )
    qt1[0,:] =(1./(9.*g))*(u_l + 2.*np.sqrt(g*h_l))**2 
    qt1[1,:] = qt1[0,:]*(u_l + 2.*(np.sqrt(g*h_l) - np.sqrt(g*qt1[0,:])))
    qt2[0,:] =(1./(9.*g))*(u_r - 2.*np.sqrt(g*h_r))**2
    qt2[1,:] = qt2[0,:]*(u_r + 2.*(np.sqrt(g*qt2[0,:]) - np.sqrt(g*h_r))) 

    # Compute q_m and associated eigenvalues
    q_m = np.zeros( (num_eqn, num_rp ) )
    q_m[0,:], q_m[1,:] = h_m, h_m*u_m
    sm[0,:] = q_m[1,:]/q_m[0,:] - np.sqrt(g*q_m[0,:])
    sm[1,:] = q_m[1,:]/q_m[0,:] + np.sqrt(g*q_m[0,:])

    # Compute waves
    wave[:,0,:] = q_m - q_l
    wave[:,1,:] = q_r - q_m

    # Evaluate q at the interface
    q = 0.5*(q_l + q_r) 
    q[:,all_shock] = q_r[:, all_shock] * (s[1,all_shock] <= 0.0) \
        + q_l[:,all_shock] * (s[0,all_shock] >= 0.0) \
        + q_m[:,all_shock] * (s[0,all_shock] < 0.0) * (0.0 < s[1,all_shock])
    q[:,one_rar] = (q_m[:,one_rar] * (sm[0,one_rar] <= 0.0) \
        + qt1[:,one_rar] * (sm[0,one_rar] >= 0.0)) * (s[0,one_rar] <= 0.0) \
        * (0.0 <= s[1,one_rar]) + q_r[:,one_rar] * (s[1,one_rar] < 0.0) \
        + q_l[:,one_rar] * (s[0,one_rar] > 0.0)
    q[:,two_rar] = (q_m[:,two_rar] * (sm[1,two_rar] >= 0.0) + qt2[:,two_rar] \
        * (sm[1,two_rar] < 0.0)) * (s[0,two_rar] <= 0.0) \
        * (0.0 <= s[1,two_rar]) + q_r[:,two_rar] * (s[1,two_rar] < 0.0) \
        + q_l[:,two_rar] * (s[0,two_rar] > 0.0)
    q[:,all_rar] = q_m[:,all_rar] * (sm[0,all_rar] <= 0.0) \
        * (0.0 <= sm[1,all_rar]) + qt1[:,all_rar] * (sm[0,all_rar] > 0.0) \
        * (s[0,all_rar] <= 0.0) + qt2[:,all_rar] * (sm[1,all_rar] < 0.0) \
        * (s[1,all_rar] >= 0.0) + q_r[:,all_rar] * (s[1,all_rar] < 0.0) \
        + q_l[:,all_rar]*(s[0,all_rar] > 0.0)

    # Compute fluctuations amdq = f(q) and apdq = -f(q)
    f = np.zeros( (num_eqn, num_rp) )
    f[0,:] = q[1,:]
    f[1,:] = ((q[1,:])**2)/q[0,:] + 0.5*g*(q[0,:])**2
    amdq, apdq  = f, -f

    return wave, s, amdq, apdq
