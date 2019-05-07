clear all
close all
clc

% Samples will be nSamples x nDimensions

% Results stored as trial(iTrial).samples(iTime, iRepetition)
% where iTrial is the current trial, which controls parameter values
% iTime is the current time step in the dynamical system
% and iRepetition is the current repetition handling the stochastic evolution

nTrials = 500; nTime = 30; nRepetitions = 500;

trial(10).samples = zeros(nTime, nRepetitions);

for iTrial = 1:nTrials
    trial(iTrial).samples = zeros(nTime, nRepetitions);
    params = hyperSamples(1);                                               % Generate parameters for this trial   
    trial(iTrial).samples(1,:) = normrnd(...
        0,(1/(1 - params(1,2)^2))^0.5, [1, nRepetitions]...                 % initialize the repetitions for each trial
        );
    for iTime = 2:nTime
        trial(iTrial).samples(iTime, :) = evolve(...
            trial(iTrial).samples(iTime - 1, :), params...                  % Perform the state evolution for each repetition simultaneously
            );
    end
end

samples = unpackWhole(trial);
[a,b] = size(samples);
samples = reshape(samples, [1, a*b]);
min(samples)
max(samples)

maxVar = var(samples);
maxStd = std(samples);
xx = linspace(min(samples), max(samples), 1000);
stdSet = maxStd*zones(size(samples));

figure(1)
histogram(samples, 'Normalization', 'pdf')
hold on
plot(xx, normpdf(xx, 0, maxStd), 'red', 'LineWidth', 3.0)
hold on
plot(xx, mixturePdf(xx, samples, maxStd*ones(size(samples))), 'blue', 'LineWidth', 3.0)


function p = mixturePdf(xx, muSet, stdSet)
    p = zeros(size(xx));
    for iDist = 1:length(muSet)
        p = p + normpdf(xx, muSet(iDist), stdSet(iDist));
    end
    p = p./iDist;
end

function samps = unpackWhole(trialObj)
    nTrials = length(trialObj);
    samps = trialObj(1).samples;
    for iTrial = 2:nTrials
        samps = [samps, trialObj(iTrial).samples];
    end
end


function inits = hyperSamples(nSamp)
    sigM = 0.01;
    mu = normrnd(0, sigM, [nSamp,1]);
    phi = betarnd(10, 1, [nSamp,1]).*2 - ones([nSamp, 1]);
    inits = [mu, phi];
end

function state = evolve(state, params)
    mu = params(1,1);
    phi = params(1,2);
    state = mu + phi*(state - mu) + normrnd(0,1, size(state));
end

function measurement = forwardModel(state)
    %measurement = normrnd(0,1, size(state)).*exp(state/2.0);
    measurement = exp(state/2.0);
end