function Nopt = optimizeMarkovChain(p,Nmin,P0,c)
% OPTIMIZEMARKOVCHAIN Compute the optimal AAD Markov chain.
%   Nopt = OPTIMIZEMARKOVCHAIN(p,Nmin,P0,c) computes the optimal number of
%   states for an AAD Markov chain based on the minimal states Nmin and
%   lower bound c of the P0-confidence interval.
%
%   Inputs:
%       p [DOUBLE]: the accuracies
%       Nmin [INTEGER]: minimal number of states
%       P0 [DOUBLE]: confidence level
%       c [DOUBLE]: lower bound confidence interval

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Input processing
[p,Nmin,P0,c] = est_utils.processInputs(p,Nmin,P0,c);
Nopt = zeros(size(p));

%% Compute optimal number of states
for l = 1:length(p)
    if p(l) == 0.5
        Nopt(l) = Inf;
    else
        N = Nmin(l);
        while stateToAmpl(lbCfdInt(p(l),N,P0(l)),N) < c(l)
            N = N+1;
        end
        Nopt(l) = N;
    end
end
end