# Takuya Ito
# 11/29/2018

# Balloon-Windkessel model for hemodynamic response modeling
# Equations and default parameters taken from: 
# Friston KJ, Harrison L, Penny W (2003) Dynamic causal modelling. Neuroimage.

import numpy as np


def balloonWindkessel(z, sampling_rate, alpha=0.32, kappa=0.65, gamma=0.41, tau=0.98, rho=0.34, V0=0.02):
    """
    Computes the Balloon-Windkessel transformed BOLD signal
    Numerical method (for integration): Runge-Kutta 2nd order method (RK2)

    z:          Measure of neuronal activity (space x time 2d array, or 1d time array)
    sampling_rate: sampling rate, or time step (in seconds)
    alpha:      Grubb's exponent
    kappa:      Rate of signal decay (in seconds)
    gamma:      Rate of flow-dependent estimation (in seconds)
    tau:        Hemodynamic transit time (in seconds)
    rho:        Resting oxygen extraction fraction
    V0:         resting blood vlume fraction

    RETURNS:
    BOLD:       The transformed BOLD signal (from neural/synaptic activity)
    s:          Vasodilatory signal
    f:          blood inflow
    v:          blood volume
    q:          deoxyhemoglobin content
    """

    if z.ndim==2:
        timepoints = z.shape[1]
    else:
        timepoints = len(z)
        z.shape = (1,len(z))

    dt = sampling_rate

    # Constants
    k1 = 7*rho
    k2 = 2
    k3 = 2*rho - 0.2

    # Create lambda function to calculate E, flow
    E = lambda x: 1.0 - (1.0 - rho)**(1.0/x) # x is f, in this case
    # Create lambda function to calculate y, the BOLD signal
    y = lambda q1,v1: V0 * (k1*(1.0-q1) + k2*(1.0 - q1/v1) + k3*(1.0 - v1))


    # initialize empty matrices to integrate through
    BOLD = np.zeros(z.shape)
    s = np.zeros(z.shape) # vasodilatory signal
    f = np.zeros(z.shape) # blood inflow
    v = np.zeros(z.shape) # blood volume
    q = np.zeros(z.shape) # deoxyhemoglobin content

    # Set initial conditions
    s[:,0] = 0.0
    f[:,0] = 1.0
    v[:,0] = 1.0
    q[:,0] = 1.0
    BOLD[:,0] = y(q[:,0], v[:,0])

    ## Obtain mean value of z, and then calculate steady state of variables prior to performing HRF modeling
    z_mean = np.mean(z,axis=1)

    # Run loop until an approximate steady state is reached 
    for t in range(timepoints-1):
    
        # 1st order increments (regular Euler)
        #s_k1 = z_mean - (1.0/kappa)*s[:,t] - (1.0/gamma)*(f[:,t] - 1.0)
        s_k1 = z_mean - (kappa)*s[:,t] - (gamma)*(f[:,t] - 1.0)
        f_k1 = s[:,t]
        v_k1 = (f[:,t] - v[:,t]**(1.0/alpha))/tau
        q_k1 = (f[:,t]*E(f[:,t])/rho - (v[:,t]**(1.0/alpha)) * q[:,t]/v[:,t])/tau

        # Compute intermediate values (Euler method)
        s_a = s[:,t] + s_k1*dt
        f_a = f[:,t] + f_k1*dt
        v_a = v[:,t] + v_k1*dt
        q_a = q[:,t] + q_k1*dt

        # 2nd order increments (RK2 method)
        #s_k2 = z_mean - (1.0/kappa)*s_a - (1.0/gamma)*(f_a - 1.0)
        s_k2 = z_mean - (kappa)*s_a - (gamma)*(f_a - 1.0)
        f_k2 = s_a
        v_k2 = (f_a - v_a**(1.0/alpha))/tau
        q_k2 = (f_a*E(f_a)/rho - (v_a**(1.0/alpha)) * q_a/v_a)/tau

        # Compute RK2 increment
        s[:,t+1] = s[:,t] + (.5*(s_k1+s_k2))*dt
        f[:,t+1] = f[:,t] + (.5*(f_k1+f_k2))*dt
        v[:,t+1] = v[:,t] + (.5*(v_k1+v_k2))*dt
        q[:,t+1] = q[:,t] + (.5*(q_k1+q_k2))*dt

        BOLD[:,t+1] = y(q[:,t+1], v[:,t+1])

        # If an approximate steady state is reached, quit.
        # We know HRF is at least 10 seconds, so make sure we wait at least 10 seconds until identifying a 'steady state'
        if (t*dt)>10 and np.sum(np.abs(BOLD[:,t+1]-BOLD[:,t]))==0: break

    ## After identifying steady state, re-initialize to run actual simulation
    s[:,0] = s[:,t+1]
    f[:,0] = f[:,t+1]
    v[:,0] = v[:,t+1]
    q[:,0] = q[:,t+1]
    BOLD[:,0] = y(q[:,t+1], v[:,t+1])

    for t in range(timepoints-1):
    
        # 1st order increments (regular Euler)
        #s_k1 = z[:,t] - (1.0/kappa)*s[:,t] - (1.0/gamma)*(f[:,t] - 1.0)
        s_k1 = z[:,t] - (kappa)*s[:,t] - (gamma)*(f[:,t] - 1.0)
        f_k1 = s[:,t]
        v_k1 = (f[:,t] - v[:,t]**(1.0/alpha))/tau
        q_k1 = (f[:,t]*E(f[:,t])/rho - (v[:,t]**(1.0/alpha)) * q[:,t]/v[:,t])/tau

        # Compute intermediate values (Euler method)
        s_a = s[:,t] + s_k1*dt
        f_a = f[:,t] + f_k1*dt
        v_a = v[:,t] + v_k1*dt
        q_a = q[:,t] + q_k1*dt

        # 2nd order increments (RK2 method)
        #s_k2 = z[:,t+1] - (1.0/kappa)*s_a - (1.0/gamma)*(f_a - 1.0)
        s_k2 = z[:,t+1] - (kappa)*s_a - (gamma)*(f_a - 1.0)
        f_k2 = s_a
        v_k2 = (f_a - v_a**(1.0/alpha))/tau
        q_k2 = (f_a*E(f_a)/rho - (v_a**(1.0/alpha)) * q_a/v_a)/tau

        # Compute RK2 increment
        s[:,t+1] = s[:,t] + (.5*(s_k1+s_k2))*dt
        f[:,t+1] = f[:,t] + (.5*(f_k1+f_k2))*dt
        v[:,t+1] = v[:,t] + (.5*(v_k1+v_k2))*dt
        q[:,t+1] = q[:,t] + (.5*(q_k1+q_k2))*dt

        BOLD[:,t+1] = y(q[:,t+1], v[:,t+1])


    return BOLD, s, f, v, q



