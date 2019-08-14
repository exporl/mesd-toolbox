%% Demonstration of the MESD toolbox to compute the MESD metric as AAD performance metric 
%
% This file demonstrates how the MESD metric can be computed based on
% evaluated performance points.
%
% NOTES: 
%   - check that the 'mesd-toolbox' directory is on the MATLAB path.
%
% The MESD is computed in four steps:
% 1. Construction of the performance curve by interpolating through the evaluated (on real EEG and audio data) (decision window length,accuracy)-points.
% 2. Optimization of the Markov chain in the number of states N for each sampled point tau on the performance curve.
% 3. Computation of the exected switch duration (ESD) per sampled tau and corresponding optimal number of states N.
% 4. The MESD is equal to the minimal ESD over all sampled decision window lengths.
%
% If this method has been useful for you, please cite the following:
% [1] S. Geirnaert, T. Francart, and A. Bertrand, “An Interpretable Performance Metric for Auditory Attention Decoding Algorithms in a Context of Neuro-Steered Gain Control,” August 2019, Internal Report
% [2] S. Geirnaert, T. Francart, and A. Bertrand, “A New Metric to Evaluate Auditory Attention Detection Performance Based on a Markov Chain,” Accepted for publication in Proc. European Signal Processing Conference (EUSIPCO), A Coruña, Spain, Sept. 2019
% [3] S. Geirnaert, T. Francart, and A. Bertrand, “MESD toolbox,” August 2019, Available online, URL: https://github.com/exporl/mesd-toolbox
%
% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

clear; close all; clc;

%% Load some dummy evaluated (decision time,accuracy)-points of three decoders

load('demoData.mat');
nbDec = 3; % number of decoders in dummy data
tau1 = tau{1}; p1 = p{1}; % the algorithm flow is first shown on one of the three decoders

%% 1. Interpolate the performance curves and display them
% The performance curve is constructed by linear interpolating through
% the evaluated points.

[tauIp,pIp] = interpolatePerfCurve(tau1,p1);

% Plot the performance curve
figure;
subplot(1,3,1);
plot(tau1,p1,'-o','linewidth',2);
xlabel('\tau [s]');
ylabel('accuracy');

%% 2. Optimize the Markov chain per sampled point on the performance curve of each decoder
% Per tau on the performance curve, a Markov chain is optimized via the
% lower bound of the P0-confidence interval. The arguments consist of the
% sampled accuracies and hyperparameter Nmin (default: 5), P0 
% (default: 0.8) and c (default: 0.65).

Nmin = 5; P0 = 0.8; c = 0.65; % the hyperparameters
Nopt = optimizeMarkovChain(pIp,Nmin,P0,c);

% Plot the optimal number of states in function of tau
subplot(1,3,2);
plot(tauIp,Nopt,'-','linewidth',2);
xlabel('\tau [s]');
ylabel('Nopt');

%% 3. Compute the expected switch duration per sampled point
% Per point on the performance curve, the expected switch duration is computed, using
% the optimized Markov chain. First, the target state is computed, whereafter 
% it can be used to compute the expected switch duration.

kc = targetState(Nopt,c);
esd = emtt(tauIp,pIp,kc);

% Plot the target state and the expected switch duration in function of tau
subplot(1,3,2); hold on;
plot(tauIp,kc,'-','linewidth',2);
ylabel('Nopt/kc');

subplot(1,3,3);
plot(tauIp,esd,'-','linewidth',2);
xlabel('\tau [s]');
ylabel('T [s]');

%% 4. Compute the MESD
% The MESD is equal to the minimal ESD over the performance curve.
[mesd,ind] = min(esd);
Nopt = Nopt(ind); tauOpt = tauIp(ind); pOpt = pIp(ind);

% Show the optimal working point
subplot(1,3,1); hold on;
plot(tauOpt,pOpt,'kd','MarkerFaceColor','k');

subplot(1,3,2); hold on;
plot(tauOpt,Nopt,'kd','MarkerFaceColor','k');
legend('Nopt','lower bound','Optimal');

subplot(1,3,3); hold on;
plot(tauOpt,mesd,'kd','MarkerFaceColor','k');

% This can all be computed with the main-function 'computeMESD':
mesd = computeMESD(tau{1},p{1});

% The hyperparameters can be changed individually by inputting for example
% computemESD(tau,p,'Nmin',10)

%% Compare now three algorithms

% Plot the performance curves
figure; hold on;
for d = 1:nbDec
    plot(tau{d},p{d},'-o','Linewidth',2);
end
xlabel('\tau [s]');
ylabel('accuracy');

% Compute the MESD
mesd = zeros(nbDec,1); Nopt = zeros(nbDec,1); tauOpt = zeros(nbDec,1); pOpt = zeros(nbDec,1);
for d = 1:nbDec
   [mesd(d),Nopt(d),tauOpt(d),pOpt(d)] = computeMESD(tau{d},p{d});
end

% Show the optimal working points
for d = 1:nbDec
    plot(tauOpt(d),pOpt(d),'kd','MarkerFaceColor','k');
end
legend('Decoder_1','Decoder_2','Decoder_3');

% Display results
mesd = table(mesd,Nopt,tauOpt,pOpt,'VariableNames',{'MESD','Nopt','tauOpt','accOpt'},'RowNames',{'Decoder_1','Decoder_2','Decoder_3'})
