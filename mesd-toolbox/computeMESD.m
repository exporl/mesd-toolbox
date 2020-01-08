function [mesd,varargout] = computeMESD(tau,p,varargin)
% COMPUTEMESD Compute the minimal expected switch duration based on
% (tau,p)-points.
%   mesd = COMPUTEMESD(tau,p) computes the minimal expected switching time
%   based on the evaluated (decision window length, accuracy)-points 
%   (tau,p), using  minimally 5 states, lower bound 0.65 and confidence 
%   level 0.8.
%
%   [mesd,Nopt,tauOpt,pOpt] = COMPUTEMESD(tau,p) also returns the
%   optimal number of states Nopt, decision window length tauOpt and 
%   accuracy pOpt.
%
%   [...] = COMPUTEMESD(tau,p,'Nmin',Nmin,'P0',P0,'c',c) uses minimal 
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
assert(all(p >= 0 & p <= 1),'p should lie within [0,1].');
if any(p <= 0.5)
   warning('accuracies below 50% (and corresponding window lengths) will be removed. Be careful when interpreting the result.');
   tau(p <= 0.5) = [];
   p(p <= 0.5) = [];
end
ip = inputParser;
addParameter(ip,'Nmin',5,@(s) assert(isnumeric(s) && all(s >= 2),'Nmin shoud be larger than or equal to 2.'));
addParameter(ip,'P0',0.8,@(s) assert(isnumeric(s) && all(s > 0 & s < 1),'P0 should lie within ]0,1[.'));
addParameter(ip,'c',0.65,@(s) assert(isnumeric(s) && all(s >= 0 & s < 1),'c should lie within [0,1[.'));
parse(ip,varargin{:})
args = ip.Results;
tau = tau(:); p = p(:);

%% Interpolate performance curve
[tau,p] = interpolatePerfCurve(tau,p);

%% Optimize Markov chain per sampled point
Nopt = optimizeMarkovChain(p,args.Nmin,args.P0,args.c);

%% Compute expected switch duration per sampled point
kc = targetState(Nopt,args.c);
esd = emtt(tau,p,kc);

%% Compute optimal ESD
[mesd,ind] = min(esd);
if tau(ind) == min(tau) || tau(ind) == max(tau)
    warning('the optimal decision length lies at the boundary of your interval. Consider taking a larger interval of decision window lengths.');
end
varargout{1} = Nopt(ind); varargout{2} = tau(ind); varargout{3} = p(ind);
end