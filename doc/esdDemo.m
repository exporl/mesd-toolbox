%% Demonstration of the ESD toolbox to compute the ESD metric as AAD performance metric 
%
% This file demonstrates how the ESD metric can be computed based on
% evaluated performance points.
%
% NOTES: 
%   - check that the 'esd-toolbox' directory is on the MATLAB path.
%
% The ESD is computed in four steps:
% 1. Construction of the performance curve by interpolating through the evaluated (on real EEG and audio data) (decision time,accuracy)-points.
% 2. Optimization of the Markov chain in the number of states N for each sampled point tau on the performance curve.
% 3. Computation of the expected Markov transit time T per sampled tau and corresponding optimal number of states N.
% 4. The ESD is equal to the minimal transit time over all evaluated transit times.
%
% If this method has been useful for you, please cite the following:
% [1] S. Geirnaert, T. Francart, and A. Bertrand, “Expected Switch Duration: an Interpretable Performance Metric to Evaluate Neural Decoders for Auditory Attention Detection,” March 2019, Internal Report
% [2] S. Geirnaert, T. Francart, and A. Bertrand, “A New Metric to Evaluate Auditory Attention Detection Performance Based on a Markov Chain,” Accepted for publication in Proc. European Signal Processing Conference (EUSIPCO), A Coruña, Spain, Sept. 2019
% [3] S. Geirnaert, T. Francart, and A. Bertrand, “ESD toolbox,” March 2019, Available online, URL: https://github.com/exporl/esd-toolbox
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

%% 3. Compute the expected Markov transit time per sampled point
% Per point on the performance curve, the expected Markov transit time is computed, using
% the optimized Markov chain. First, the lower bound of the P0-confidence
% interval is computed, whereafter it can be used to compute the expected Markov transit
% time.

k = ceil(c.*(Nopt-1)+1);
T = emtt(tauIp,pIp,k);

% Plot the lower bound of the P0-confidence interval and the expected 
% Markov transit time in function of tau
subplot(1,3,2); hold on;
plot(tauIp,k,'-','linewidth',2);
ylabel('Nopt/k');

subplot(1,3,3);
plot(tauIp,T,'-','linewidth',2);
xlabel('\tau [s]');
ylabel('T [s]');

%% 4. Compute the ESD
% The ESD is equal to the minimal expected Markov transit time over the performance curve.
[esd,ind] = min(T);
Nopt = Nopt(ind); tauOpt = tauIp(ind); pOpt = pIp(ind);

% Show the optimal working point
subplot(1,3,1); hold on;
plot(tauOpt,pOpt,'kd','MarkerFaceColor','k');

subplot(1,3,2); hold on;
plot(tauOpt,Nopt,'kd','MarkerFaceColor','k');
legend('Nopt','lower bound','Optimal');

subplot(1,3,3); hold on;
plot(tauOpt,esd,'kd','MarkerFaceColor','k');

% This can all be computed with the main-function 'ESD':
esd = computeESD(tau{1},p{1});

% The hyperparameters can be changed individually by inputting for example
% computeESD(tau,p,'Nmin',10)

%% Compare now three algorithms

% Plot the performance curves
figure; hold on;
for d = 1:nbDec
    plot(tau{d},p{d},'-o','Linewidth',2);
end
xlabel('\tau [s]');
ylabel('accuracy');

% Compute the ESD
esd = zeros(nbDec,1); Nopt = zeros(nbDec,1); tauOpt = zeros(nbDec,1); pOpt = zeros(nbDec,1);
for d = 1:nbDec
   [esd(d),Nopt(d),tauOpt(d),pOpt(d)] = computeESD(tau{d},p{d});
end

% Show the optimal working points
for d = 1:nbDec
    plot(tauOpt(d),pOpt(d),'kd','MarkerFaceColor','k');
end
legend('Decoder_1','Decoder_2','Decoder_3');

% Display results
T = table(esd,Nopt,tauOpt,pOpt,'VariableNames',{'ESD','Nopt','tauOpt','accOpt'},'RowNames',{'Decoder_1','Decoder_2','Decoder_3'})
