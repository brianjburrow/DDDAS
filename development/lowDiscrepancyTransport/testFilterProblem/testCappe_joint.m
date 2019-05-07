% generate a uniform proposal set over the problem space
clear all
clc
close all
%% Params to move to function definition
numSamples = 1000;

%% Params
nRepeats = 1000;
nTime = 70;
%% Generate joint samples

counter = 1;
for iRepeat = 1:nRepeats
    currentState = normrnd(0, 10^0.5);
    for iTime = 1:nTime
        states(counter) = cappe_benchmark.StateTransition(currentState, iTime);
        currentState = states(counter);
        counter = counter + 1;
    end
end


%% Set Bounds

lowBound = min(states);
upBound = max(states);

%% generate sobol set
propSamples = states;

output = cappe_benchmark.GetDataC(propSamples);



subplot(2,1,1)
scatter(propSamples, output)
subplot(2,1,2)
histogram(propSamples)