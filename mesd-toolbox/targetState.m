function kc = targetState(N,c)
% LBCFDINT Compute the target state of the AAD Markov chain.
%   kc = TARGETSTATE(N,c) computes the target state of the AAD Markov chain
%   as the first state corresponding to a relative gain x > c, with comfort
%   level c.
%
%   Inputs:
%       N [INTEGER]: number of states
%       c [DOUBLE]: comfort level (default: 0.65)

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and processing inputs
[N,c] = esd_utils.processInputs(N,c);

%% Compute target state
kc = ceil(c.*(N-1)+1);
end