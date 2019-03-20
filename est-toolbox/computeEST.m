function [est,varargout] = computeEST(tau,p,varargin)
% COMPUTEEST Compute the expected switching time based on (tau,p)-points.
%   est = COMPUTEEST(tau,p) computes the expected swithcing time based on the 
%   evaluated (decision time, accuracy)-points (tau,p), using minimally 5 
%   states, lower bound 0.65 and confidence level 0.9.
%
%   [est,Nopt,tauOpt,pOpt] = COMPUTEEST(tau,p) also returns the
%   optimal number of states Nopt, decision time tauOpt and accuracy
%   pOpt.
%
%   [...] = EST(tau,p,'Nmin',Nmin,'P0',P0,'c',c) uses minimal number of
%   states Nmin, confidence level P0 and lower bound c.
%
%   Inputs:
%       tau [DOUBLE]: evaluated decision lengths
%       p [DOUBLE]: evaluated accuracies
%       (optional) 'Nmin', Nmin [INTEGER]: minimal number of states (default: 5)
%       (optional) 'P0', P0 [DOUBLE]: confidence level (default: 0.9)
%       (optional) 'c', c [DOUBLE]: lower bound confidence interval (default: 0.65)

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% Asserts and input processing
[tau,p] = est_utils.processInputs(tau,p);
assert(all(tau > 0),'tau should be positive.');
assert(all(p > 0.5 & p <= 1),'p should lie within ]0.5,1].');
ip = inputParser;
addParameter(ip,'Nmin',5,@(s) assert(isinteger(s) && all(s >= 2),'Nmin shoud be larger than or equal to 2.'));
addParameter(ip,'P0',0.9,@(s) assert(isnumeric(s) && all(s > 0 & s < 1),'P0 should lie within ]0,1[.'));
addParameter(ip,'c',0.65,@(s) assert(isnumeric(s) && all(s >= 0 & s < 1),'c should lie within [0,1[.'));
parse(ip,varargin{:})
args = ip.Results;
tau = tau(:); p = p(:);

%% Interpolate performance curve
[tau,p] = interpolatePerfCurve(tau,p);

%% Optimize Markov chain per sampled point
Nopt = optimizeMarkovChain(p,args.Nmin,args.P0,args.c);

%% Compute transit time per sampled point
k = lbCfdInt(p,Nopt,args.P0);
T = transitTime(tau,p,k);

%% Compute optimal MTT
[est,ind] = min(T);
varargout{1} = Nopt(ind); varargout{2} = tau(ind); varargout{3} = p(ind);
end