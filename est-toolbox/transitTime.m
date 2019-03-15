function T = transitTime(tau,p,k)
% TRANSITTIME Compute the expected transit time to state k.
%   T = TRANSITTIME(tau,p,k) computes the transit time to state k, averaged
%   over all initial states (1,...,k-1), with given transition probability 
%   p and decision time tau.
%
%   Inputs:
%       tau [DOUBLE]: the decision time
%       p [DOUBLE]: the transition probability (0.5 < p <= 1)
%       k [INTEGER]: the target state

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and processing inputs
[tau,p,k] = est_utils.processInputs(tau,p,k);

%% Compute transit time
T = zeros(length(tau),1);
for l = 1:length(tau)
    if k(l) == 1
        T(l) = 0;
    else
        ss = flipud(ssDistribution(p(l),k(l)+1)); % note: the number of states is irrelevant, because of the normalization (conditioning)
        ss = ss./sum(ss(1:k(l)-1)); % normalization
        for i = 1:k(l)-1
            T(l) = T(l) + ss(i)*meanHittingTime(p(l),i,k(l)); % per initial state contribution
        end
        T(l) = T(l).*tau(l);
    end
end
end