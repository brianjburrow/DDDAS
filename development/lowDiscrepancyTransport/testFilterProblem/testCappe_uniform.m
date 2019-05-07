%% testCappe_uniform 
% generate a uniform proposal set over the problem space
clear all
clc
close all
%% Params to move to function definition
numSamples = 70000;

%% Params
nRepeats = 1000;
nTime = 70;
%% Generate Joint Samples

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

ss = sobolset(1);
ss = scramble(ss, 'MatousekAffineOwen');
propSamples = ss(1:numSamples,:)*(upBound - lowBound) + lowBound;

output = cappe_benchmark.GetDataC(propSamples);


%% Set simulation parameters
processNoise = 10;
state(1) = 0;
measurementNoise = 1;

%% Initialize L2O_Filter

L2O.samples = propSamples';
L2O.data = cappe_benchmark.GetDataC(0);
L2O.stateEvolution = @(xx, yy) cappe_benchmark.StateTransition(xx, yy); %xx is the state, yy is time;
L2O.oper = 0;
L2O.output = output';
L2O.dataFunc = @(xx, yy) cappe_benchmark.GetDataC(xx);
L2O.R = measurementNoise;
L2O.estimate = 0;
L2O.err = processNoise;

for idx = 2:70
    disp(idx)
    % Evolve State and collect data
    state(idx) = cappe_benchmark.StateTransition(state(idx - 1), idx);
    data = cappe_benchmark.GetData(state(idx));
    
    %Evaluate L2O_filter
    L2O.oper = idx;
    L2O.data = data;
    tic
    L2O = L2oFilter(L2O, 0);
    l2oTime(idx) = toc;
    est(idx) = L2O.estimate;
    err(idx) = L2O.err;     
end

figure(1)
plot(state, 'black', 'LineWidth', 2.0, 'DisplayName', 'Truth')
hold on
plot(est, 'blue', 'LineWidth', 2.0, 'DisplayName', 'Estimate')
legend()
xlabel("Time")
ylabel("State")
figure(2)
plot(state - est, 'red', 'LineWidth', 2.0, 'DisplayName', 'Error')
xlabel("Time")
ylabel("Error")
legend()