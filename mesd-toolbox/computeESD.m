function [esd,varargout] = computeESD(tau,p,varargin)
% COMPUTEESD Compute the expected switch duration for given (tau,p)-points.
%   esd = COMPUTEESD(tau,p) computes the expected switch duration
%   for given (decision window length, accuracy)-point(s) 
%   (tau,p), using  minimally 5 states, lower bound 0.65 and confidence 
%   level 0.8.
%
%   [esd,Nopt,kc] = COMPUTEESD(tau,p) also returns the
%   optimal number of states Nopt and target state kc.
%
%   [...] = COMPUTEESD(tau,p,'Nmin',Nmin,'P0',P0,'c',c) uses minimal 
%   number of states Nmin, confidence level P0 and lower bound c.
%
%   Inputs:
%       tau [DOUBLE]: evaluated decision window lengths
%       p [DOUBLE]: evaluated accuracies
%       (optional) 'Nmin', Nmin [INTEGER]: minimal number of states (default: 5)
%       (optional) 'P0', P0 [DOUBLE]: confidence level (default: 0.8)
%       (optional) 'c', c [DOUBLE]: comfort level (default: 0.65)

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and input processing
[tau,p] = esd_utils.processInputs(tau,p);
assert(all(tau > 0),'tau should be positive.');
assert(all(p > 0.5 & p <= 1),'p should lie within ]0.5,1].');
ip = inputParser;
addParameter(ip,'Nmin',5,@(s) assert(isnumeric(s) && all(s >= 2),'Nmin shoud be larger than or equal to 2.'));
addParameter(ip,'P0',0.8,@(s) assert(isnumeric(s) && all(s > 0 & s < 1),'P0 should lie within ]0,1[.'));
addParameter(ip,'c',0.65,@(s) assert(isnumeric(s) && all(s >= 0 & s < 1),'c should lie within [0,1[.'));
parse(ip,varargin{:})
args = ip.Results;
tau = tau(:); p = p(:);

%% Optimize Markov chain
Nopt = optimizeMarkovChain(p,args.Nmin,args.P0,args.c);

%% Compute expected switch duration
kc = targetState(Nopt,args.c);
esd = emtt(tau,p,kc);
varargout{1} = Nopt; varargout{2} = kc;
end