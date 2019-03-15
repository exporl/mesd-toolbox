function ss = ssDistribution(p,N)
% SSDISTRIBUTION Compute the steady-state distribution of the AAD Markov
% chain.
%   ss = SSDISTRIBUTION(p,N) computes the steady-state distribution of the
%   AAD Markov chain with transition probability p and N states.
%
%   Inputs:
%       p [DOUBLE SCALAR]: the transition probability (0 <= p <= 1)
%       N [INTEGER SCALAR]: the number of states

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

q = 1-p;
r = p./q;
ss = (r-1)/(r^N-1)*r.^((1:N)'-1);
end