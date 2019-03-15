function MHT = meanHittingTime(p,i,k)
% MEANHITTINGTIME Compute mean hitting time from state i to k.
%   MHT = MEANHITTINGTIME(p,i,k) computes the mean hitting time from
%   state i to given state k, given transition probability p ..
%
%   Inputs:
%       p [DOUBLE]: the transition probability (0 <= p <= 1)
%       i [INTEGER]: the initial state
%       k [INTEGER]: the target state

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Process inputs and assert i <= k
[p,i,k] = est_utils.processInputs(p,i,k);
assert(all(i <= k),'i should be smaller or equal to k. If you want i > k, change p -> 1-p, i -> N-i+1 and k -> N-k+1');

%% Compute MHT
q = 1-p;
r = p./q;
MHT = (k-i)./(p-q)+(p.*(r.^(-k)-r.^(-i)))./(p-q).^2;

%% Correct for boundary effects
MHT(p==0.5) = (k(p==0.5)-i(p==0.5)).*(k(p==0.5)+i(p==0.5)-1);
MHT(p==0) = inf;
end