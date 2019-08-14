% EST TOOLBOX
% Version 1.2, 14-08-2019
%
% MAIN FUNCTION
%   computeMESD           - Compute the minimal expected switch duration 
%                           based on (tau,p)-points.
%
% OPTIMIZATION MARKOV CHAIN MODEL
%   optimizeMarkovChain  - Compute the optimal AAD Markov chain.
%
% BASIC FEATURES MARKOV CHAIN MODEL
%   ssDistribution       - Compute the steady-state distribution of the AAD
%                          Markov chain.
%   stateToAmpl          - Convert a state index to a relative amplification level.
%   lbCfdInt             - Compute the lower bound of the P0-confidence
%                          interval of the AAD Markov chain.
%   evState              - Compute the expected value of the AAD Markov chain.
%   targetState          - Compute the target state of the AAD Markov chain.
%
% FUNDAMENTAL METRICS
%   meanHittingTime      - Compute the mean hitting time from state i to k.
%   emtt                 - Compute the expected Markov transit time to state k.
%                          The EMTT to state kc is equal to the ESD.
%
% CONSTRUCTION PERFORMANCE CURVE
%   interpolatePerfCurve - Interpolate the performance curve through
%                          evaluated performance points.