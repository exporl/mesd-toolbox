function x = stateToAmpl(state,N)
% STATETOAMPL Convert a state index to a relative amplification level.
%   x = STATETOAMPL convert the state index state in a AAD Markov chain
%   with N states to relative amplification level.
%
%   Inputs:
%       state [INTEGER]: the state
%       N [INTEGER]: the number of states

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and input processing
[state,N] = est_utils.processInputs(state,N);

%% Transform
x = (state-1)./(N-1);
end