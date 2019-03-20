% EST TOOLBOX
% Version 1.0, 06-03-2019
%
% MAIN FUNCTION
%   computeEST           - Compute the expected switching time based on (tau,p)-points.
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
%
% FUNDAMENTAL METRICS
%   meanHittingTime      - Compute the mean hitting time from state i to k.
%   transitTime          - Compute the expected transit time to state k.
%
% CONSTRUCTION PERFORMANCE CURVE
%   interpolatePerfCurve - Interpolate the performance curve through
%                          evaluated performance points.