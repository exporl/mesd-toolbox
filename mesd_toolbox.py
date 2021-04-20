# MESD TOOLBOX (Python)
# Version 1.0, 19-04-2021
#
# Python version, based on Simon Geirnaert's MESD toolbox for Matlab.
#
# Author: Debora Fieberg, KU Leuven
# Department of Neurosciences, Research Group Experimental Oto-rhino-laryngology &
# Department of Electrical Engineering (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data Analytics
# Correspondence: debora.fieberg@kuleuven.be
#
#
# MAIN FUNCTIONS
#   compute_ESD                     - Compute the expected switch duration for given (tau,p)-point(s).
#   compute_MESD                    - Compute the minimal expected switch duration based on (tau,p)-points.
#   optimize_Markov_chain           - Compute the optimal AAD Markov chain.
#   compute_targetState_and_EMTT    - Compute the target state and the respective expected Markov transit time.
#   interpolate_performance_curve   - Interpolate the performance curve through evaluated performance points.


import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import warnings


def compute_ESD(tau, p, N_min=5, P0=0.8, c=0.65):
    ''' Compute the expected switch duration (ESD) for given (decision window length, accuracy)-points (tau,p), 
        by default using  minimally 5 states, confidence level 0.8, and lower bound 0.65.

        Usage:
            ESD, N_opt, k_c = compute_ESD(tau,p,N_min,P0,c)
            ESD, N_opt, k_c = compute_ESD(tau,p)            # using default values for P0, N_min, and c
            ESD, *_  = compute_ESD(tau,p,N_min,c)           # use when just interested in the ESD

        Inputs:
            tau [list or np.ndarray of numeric values]: evaluated decision window lengths
            p [list or np.ndarray of numeric values]: evaluated accuracies
            (optional) N_min [int]: minimal number of states (default: 5)
            (optional) P0 [float]: confidence level (default: 0.8)
            (optional) c [float]: comfort level (default: 0.65)

        Outputs:
            ESD : the expected switch durations
            N_opt : the optimal numbers of states in the Markov chain
            k_c : target states
    '''

    # Asserts and input processing
    tau, p, N_min, P0, c = check_and_process_inputs(tau, p, N_min, P0, c)

    # Optimize Markov chain
    N_opt = optimize_Markov_chain(tau, p, N_min, P0, c)

    # Compute expected switch duration
    ESD, k_c = compute_targetState_and_EMTT(tau, p, N_opt, c)
    
    return ESD, N_opt, k_c




def compute_MESD(tau, p, N_min=5, P0=0.8, c=0.65):
    ''' Compute the minimal expected switch duration based on (tau,p)-points.
    mesd = compute_MESD(tau,p) computes the minimal expected switching time
    based on the evaluated (decision window length, accuracy)-points (tau,p), 
    using  minimally 5 states, lower bound 0.65 and confidence level 0.8.

    Usage:
        mesd, N_opt, tau_opt, p_opt = compute_MESD(tau,p,N_min,P0,c)
        mesd, N_opt, tau_opt, p_opt = compute_MESD(tau,p)            # using default values for P0, N_min, and c
        mesd, *_  = compute_MESD(tau,p,N_min,P0,c)                   # use when just interested in the mesd

    Inputs:
        tau [list or np.ndarray of numeric values]: evaluated decision window lengths
        p [list or np.ndarray of numeric values]: evaluated accuracies
        (optional) N_min [INTEGER]: minimal number of states (default: 5)
        (optional) P0 [DOUBLE]: confidence level (default: 0.8)
        (optional) c [DOUBLE]: comfort level (default: 0.65)

    Outputs:
        mesd : the minimal expected switch duration
        N_opt : the optimal number of states
        tau_opt : decision window length optimal w.r.t. the MESD
        p_opt : accuracy optimal w.r.t. the MESD
    '''
    
    # Asserts and input processing
    tau, p, N_min, P0, c = check_and_process_inputs(tau, p, N_min, P0, c)
    
    # Interpolate performance curve
    tau, p = interpolate_performance_curve(tau,p)

    # Optimize Markov chain per sampled point
    N_opt = optimize_Markov_chain(tau, p, N_min, P0, c)
    
    # Compute expected switch duration per sampled point
    ESDs, kc = compute_targetState_and_EMTT(tau, p, N_opt, c)

    # Compute optimal ESD
    mesd = np.min(ESDs); index = np.argmin(ESDs)
    if tau[index] == np.min(tau) or tau[index] == np.max(tau):
        warnings.warn("The optimal decision length lies at the boundary of your interval. Consider taking a larger interval of decision window lengths.")

    return mesd, N_opt[index], tau[index], p[index]




def optimize_Markov_chain(tau, p, N_min=5, P0=0.8, c=0.65):
    '''Compute the optimal AAD Markov chain.
        N_opt = optimize_Markov_chain(p,N_min,P0,c) computes the optimal number of states for an AAD Markov chain based 
        on the minimal states Nmin and lower bound c of the P0-confidence interval.
        
        Usage:
            N_opt = optimize_Markov_chain(tau, p, N_min, P0, c)
            N_opt = optimize_Markov_chain(tau, p)                        # using default values for P0, N_min, and c

        Inputs:
            tau [list or np.ndarray of numeric values]: evaluated decision window lengths
            p [list or np.ndarray of numeric values]: evaluated accuracies
            N_min [int]: minimal number of states
            P0 [float]: confidence level
            c [float or int]: lower bound confidence interval

        Outputs:
            N_opt: the optimal number of states for an AAD Markov chain 
    '''

    # Input processing
    #tau, p ,N_min, P0, c = check_and_process_inputs(tau, p, N_min, P0, c)

    # Compute optimal number of states
    N_opt = []
    for index, wl in enumerate(tau):
        acc = p[index]; r = acc/(1-acc)
        if acc == 0.5:
            N = np.inf
        else:
            N = N_min; kbar = np.floor(np.log( r**N *(1-P0) +P0 )/np.log(r) +1)
            while (kbar-1)/(N-1) < c:
                N+=1; kbar = np.floor(np.log( r**N *(1-P0) + P0 )/np.log(r) +1)
                # TODO: correct for boundary effects
        N_opt.append(N)

    return N_opt




def compute_targetState_and_EMTT(tau, p, Nopt, c = 0.65):
    ''' Compute the target state kc and the expected Markov transit time (EMTT).
        T, kc = compute_targetState_and_EMTT(tau, p, Nopt, P0, c) computes the target state kc
        and the EMTT to state kc, averaged over all initial states (1,...,k-1), with given transition
        probability p and decision time tau. The EMTT with target state kc is equal to the ESD.

        Usage:
            ESDs, kc = compute_targetState_and_EMTT((tau, p, Nopt, c)

        Inputs:
            tau [list]: the decision time
            p [list]: the transition probability (0.5 < p <= 1)
            Nopt [list]: the optimal number of states for the AAD Markov chains 
            (optional) c [float]: comfort level (default: 0.65)

        Outputs:
            ESDs : list of expected switch durations
            kc : list of target states
    '''

    assert len(tau) == len(p), "Input size mismatch"
    assert len(tau) == len(Nopt), "Input size mismatch"

    ESDs = []; kc = []

    for index, wl in enumerate(tau):

        acc = p[index]; r = acc/(1-acc)

        # Compute the target state
        k = int(np.ceil(c*(Nopt[index]-1)+1)) 
        kc.append(k) 
        
        # Compute the expected Markov transit time (EMTT) to state k
        if k == 1:
            ESD = 0
        else:
            summand = 0
            for i in range(1,k):
                h = (k-i)/(2*acc-1) + (acc* (r**(-k)-r**(-i)) )/((2*acc-1)**2) # mean hitting time
                summand += r**(-i) * h
            ESD = wl * (r**(k+1)-r**k)/(r**k-r) * summand
            ESDs.append(ESD)

    return ESDs, kc




def interpolate_performance_curve(tau,p):
    ''' Linearly interpolate the performance curve through the evaluated points.
          [tau,p] = interpolate_performance_curve(tau,p) linearly interpolates the performance through 
          evaluated (decision window length, accuracy)-points (tau,p).
    
    Usage:
        tau_Ip, p_Ip = interpolate_performance_curve(tau,p)

    Inputs:
        tau [list of floats]: evaluated decision window lengths
        p [list of floats]: evaluated accuracies

    Outputs:
        tau_Ip : interpolated window lengths
        p_Ip : interpolated accuracies
    '''

    assert len(tau) == len(p), "Input size mismatch"

    #  Sort lists according to window lengths tau:
    zipped_lists = zip(tau, p)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    tau, p = [ list(tuple) for tuple in  tuples]

    # Linearly interpolate performance curve through evaluated points:
    f = interpolate.interp1d(tau, p)
    tau_Ip = np.linspace(min(tau),max(tau),1000)
    p_Ip = f(tau_Ip)

    return tau_Ip, p_Ip





def process_input(x):
    '''Check that and input x is of the correct type and transform it to numpy array.'''

    # input should be a list or a 1D array:
    assert isinstance(x,(list,np.ndarray)), "Window lengths and accuracies have to be given as lists or 1D numpy arrays."
    if x is np.ndarray:
        assert len(np.shape(x)) == 1, "Window lengths and accuracies have to be given as lists or 1D numpy arrays."
    
    # all entries should be numeric:
    #assert all( isinstance(ele, (int, float)) for ele in x), "Window lengths and accuracies should be numeric values (int or float)" 
    return np.array(x)




def check_and_process_inputs(tau, p, N_min, P0, c):
    
    tau = process_input(tau); p = process_input(p)
    assert len(tau) == len(p), "Input size mismatch"
    
    assert np.all(tau > 0), "tau should be positive."
    assert np.all(p >= 0) and np.all(p<=1), "p should lie within [0,1]."
    assert np.any(p > 0.5), "All accuracies are below 50%. The MESD cannot be computed."
    if np.any(p <= 0.5):
        warnings.warn("Accuracies below 50% (and corresponding window lengths) will be removed. Be careful when interpreting the result.")
        tau = tau[np.where(p>0.5)]; p = p[np.where(p>0.5)]
    
    assert type(N_min) is int, "N_min shoud be an integer."
    assert N_min >= 2, "N_min shoud be larger than or equal to 2."

    assert isinstance(P0, float), "P0 should be a float value within ]0,1[."
    assert P0>0 and P0<1, "P0 should be a value within ]0,1[."

    assert isinstance(c, (float, int)), "c should be a value within [0,1[."
    assert c>=0 and c<1, "c should be a value within [0,1[."

    return list(tau), list(p), N_min, P0, c





