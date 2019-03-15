function [tau,p] = interpolatePerfCurve(tau,p)
% INTERPOLATEPERFCURVE Interpolate the performance curve through evaluated
% points.
%   [tau,p] = INTERPOLATEPERFCURVE(tau,p) linearly interpolates the
%   performance through evaluated (decision time, accuracy)-points (tau,p).
%
%   Inputs:
%       tau [DOUBLE]: evaluated decision lengths
%       p [DOUBLE]: evaluated accuracies

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and input processing
assert(length(tau)==length(p),'Input size mismatch');
tau = tau(:); p = p(:);

%% Compute performance curve
tauTemp = linspace(min(tau),max(tau),1e3)';
[tau,S] = sort(tau,'ascend');
p = p(S);
p = interp1q(tau,p,tauTemp);
tau = tauTemp;
end