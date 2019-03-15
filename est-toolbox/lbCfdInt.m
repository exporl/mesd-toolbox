function k = lbCfdInt(p,N,P0)
% LBCFDINT Compute the lower bound of the P0-confidence interval of the AAD
% Markov chain.
%   k = LBCFDINT(p,N,P0) computes the lower bound of the P0-confidence
%   interval of the AAD Markov chain with transition probability p and N
%   states.
%
%   Inputs:
%       p [DOUBLE]: accuracy
%       N [INTEGER]: number of states
%       P0 [DOUBLE]: confidence level

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and processing inputs
[p,N,P0] = est_utils.processInputs(p,N,P0);

%% Compute lower bound
r = p./(1-p);
k = floor(log(P0+(1-P0).*r.^N)./log(r)+1);

%% Correct for boundary effects
k(p==1) = N(p==1);
k(p==0.5) = floor(N(p==0.5).*(1-P0(p==0.5))+1);
end