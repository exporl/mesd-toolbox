function E = evState(p,N)
% EVSTATE Compute the expected value of the AAD Markov chain.
%   E = EVSTATE(p,N) computes the expected value of the
%   AAD Markov chain with N states and transition probability p.
%
%   Input parameters:
%       p [DOUBLE]: the transition probability (probability of success)
%       N [INTEGER]: the number of states in the Markov chain

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Process inputs
[p,N] = est_utils.processInputs(p,N);

%% Compute expected value
r = p./(1-p);
E = round((N.*r.^N)./(r.^N-1)-1./(r-1));

switch p
    case 0.5
        E = (N+1)/2;
    case 1
        E = N;
end
end