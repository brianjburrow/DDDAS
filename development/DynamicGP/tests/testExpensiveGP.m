%% testDynamicGP_filteringCost.m
% Author: Brian Burrows
% This method is to test the cost, and accuracy of various surrogates
% using a fixed training set size

clear all
close all
clc


%% Problem Initialization
nUpdates = 5;
nClusters = 20;


lowerBound = 0.5;
lb = lowerBound;
upperBound = 2.5;
ub = upperBound;
librarySize = 200;
nMonteCarlo = 100;
Tmax = 70;
libraryIn = linspace(lowerBound, upperBound, librarySize)';                  % Fix the
libraryOut = GetDataC(libraryIn);

%% Filtering Initialization for Cappe Test Function
initialState = 1.5;
initialUncertainty = 0.5;
processNoise = 0.001;
measurementNoise = 0.1;
initialMeasurement = GetData(initialState);

%% Optimize GP Hyperparameters
% Initial GP parameters
ell = 0.05;                                                                 % Characteristic length scale of covariance function
sf = 3;                                                                     
sn = 0.2; 
% Set GP mean, covariance, and likelihood functions
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0; 0];
covfunc = {@covMaterniso, 3}; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.2; hyp.lik = log(sn);

% Perform optimization
noisyOutputs = GetData(libraryIn')';
hyperParam = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                 covfunc, likfunc, libraryIn, noisyOutputs); 

%% Set up instance of the dynamicGP class
% set up GPML function handle
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, ...
    meanfunc, covfunc,...
    likfunc, xx,...
    yy, zz);
% Set up initial pdf over our input space
inputSampler = @(numSamp) normpdf(initialState, initialUncertainty);
trainIn =  libraryIn(1:30:end);                                 % only place the bounds of the problem in the training set
trainOut = libraryOut(1:30:end);                                % only place the outputs corresponding to the bounds of the training set

[maxTrainingSize, ~] = size(trainIn);

dynGP1 = dynamicGP(gpHandle, trainIn,...
            trainOut, inputSampler,...
            libraryIn, libraryOut);

%% Train the initial GP
disp(maxTrainingSize - length(trainIn))
dynGP1 = dynGP1.addNewPoints(nMonteCarlo, maxTrainingSize - length(trainIn));
[muVec, stdVec] = dynGP1.predict(libraryIn);

figure(10)
niceErrorbar(libraryIn, muVec, 3*stdVec)
hold on
scatter(dynGP1.trainInput, dynGP1.trainOutput)
drawnow
figure(1)
dynGP2 = dynGP1.copy();
dynGP3 = dynGP2.copy();

dynGP1.trainInput = libraryIn(1:5:end);
dynGP1.trainOutput = libraryOut(1:5:end);


%% Initiazile Bootstrap Particle Filter
numParticles = 5000;                                                        % Set the number of particles used at each iteration.  Increase for better accuracy but higher computational cost
pdf_g = @(dat,xx) normpdf(dat, 0.05 * xx^2, 1);                             % Set up the likelihood model

PF1.resampling = 1;                                                          % Set resampling = true, this will change dynamically later
PF1.likelihood = pdf_g; 
PF1.proposalSampler.random = @(xx, dat, Tim) normrnd(...                     % Create a random number generator to sample from the state evolution function
    StateTransitionC(xx, Tim),...
    processNoise^0.5);
PF1.particles = zeros(length(initialState), numParticles);                   
PF1.weights = zeros(1, numParticles);
for idx = 1:numParticles
    data = GetData(initialState);                           %
    PF1.particles(:,idx) = PF1.proposalSampler.random(initialState, ...       % generate an initial particle
        data, 0);
    PF1.weights(idx) = PF1.likelihood(data, initialState);                    % Compute initial weights
end
PF1.weights = PF1.weights/sum(PF1.weights);                                    % Normalize the weights


PF2 = PF1;
PF3 = PF1;

states = initialState;
errs = initialUncertainty^0.5;
measurements = initialMeasurement;
for iTime = 2:Tmax
    fprintf("Processing Time Step %d \n", iTime)
    %% Perform Particle filtering
    PF1.likelihood = @(data, possibleStates) dynGP1.likelihood(data, possibleStates')';
    PF2.likelihood = @(data, possibleStates) dynGP2.likelihood(data, possibleStates')';
    PF3.likelihood = @(data, possibleStates) dynGP3.likelihood(data, possibleStates')';
    
    %% Plotting
    figure(10)
    subplot(2,2, 1)
    [muVec, stdVec] = dynGP1.predict(libraryIn);
    niceErrorbar(libraryIn, muVec, 3*stdVec)
    hold on
    ylim([-3, 6])
    xlim([lb, ub])
    xlabel('X')
    ylabel('f(X)')
    scatter(dynGP1.trainInput, dynGP1.trainOutput, 10, 'filled', 'b')
    drawnow
    hold off
    
    
    subplot(2,2, 3)
    [muVec, stdVec] = dynGP3.predict(libraryIn);
    niceErrorbar(libraryIn, muVec, 3*stdVec)
    hold on
    ylim([-3, 6])
    xlim([lb, ub])
    xlabel('X')
    ylabel('f(X)')
    scatter(dynGP3.trainInput, dynGP3.trainOutput, 10, 'filled', 'b')
    drawnow
    hold off
    % Evaluate the Bootstrap Particle Filter
    PF1.data = measurements(iTime - 1);
    PF1.oper = iTime - 1;
    tic
    PF1 = BootstrapFilter(PF1);
    a = toc;
    pf1Time(iTime - 1) = a;
    
    pf1Est(iTime - 1) = PF1.estimate;
    pf1Err(iTime - 1) = PF1.err^0.5;
    fprintf("dynGP1 took %f seconds \n", a)
    
    % Evaluate the Bootstrap Particle Filter
    PF3.data = measurements(iTime - 1);
    PF3.oper = iTime - 1;
    tic
    PF3 = BootstrapFilter(PF3);
    %% Evolve the Gaussian Process Model to be better at the next time step
    dynGP3.inputSampler = @(numSampl) StateTransition(...
        PF3.particles(:, 1:numSampl), iTime)';
    dynGP3.updateGaussianProcess_kMeans(numParticles, nUpdates, nClusters);
    updateTime3(iTime) = toc;
    
    a = toc;
    fprintf("dynGP3 took %f seconds \n", a)
    pf3Time(iTime - 1) = a;
    
    pf3Est(iTime - 1) = PF3.estimate;
    pf3Err(iTime - 1) = PF3.err^0.5;


    if iTime > 1
        subplot(2, 2, 2)
        niceErrorbar([1:iTime - 1]', zeros([iTime - 1, 1]), pf1Err(1:iTime - 1)');
        hold on
        scatter(1:iTime - 1, states(1:iTime - 1) - pf1Est(1:iTime - 1), 10, 'filled', 'r')
        ylim([-2, 2])
        xlim([0, 70])
        xlabel('Time')
        ylabel('Error')
        hold off

        subplot(2, 2, 4)
        niceErrorbar([1:iTime - 1]', zeros([iTime - 1, 1]), pf3Err(1:iTime - 1)');
        hold on
        scatter(1:iTime - 1, states(1:iTime - 1) - pf3Est(1:iTime - 1), 10, 'filled', 'r')
        ylim([-2, 2])
        xlim([0, 70])
        xlabel('Time')
        ylabel('Error')
        hold off
    end
    
    %% Evolve state and collect measurement
    if iTime < Tmax
        states(iTime) = StateTransition(states(iTime -1), ...
            iTime);
        measurements(iTime) = GetData(states(iTime));
    end
   
end


e1 = mean((states - pf1Est).^2)^0.5;
e3 = mean((states - pf3Est).^2)^0.5;

figure(1)
subplot(2,1,1)
niceErrorbar([1:Tmax-1]', zeros([Tmax-1,1]), pf1Err')
hold on
scatter([1:Tmax-1], states - pf1Est, 10, 'filled', 'k');
plot([1:Tmax-1], e1*ones([1, Tmax-1]))

subplot(2,1,2)
niceErrorbar([1:Tmax-1]', zeros([Tmax-1,1]), pf3Err')
hold on
scatter([1:Tmax-1], states - pf3Est, 10, 'filled', 'k');
plot([1:Tmax-1], e3*ones([1, Tmax-1]))

function state = StateTransitionC(state, time)
     % rescale from [0.5, 2.5] to [-60, 60]
    minSet = 0.5 * ones(size(state));
    maxSet = 2.5 * ones(size(state));
    state = 120*(state - minSet) ./ (maxSet - minSet) - 60*ones(size(state));
    state = cappe_benchmark.StateTransitionC(state, time);
    % rescale from [-60, 60] to [0.5, 2.5]
    minSet2 = -60 * ones(size(state));
    maxSet2 = 60 * ones(size(state));
    state = (2.5 - 0.5) * (state - minSet2) ./ (maxSet2 - minSet2) + 0.5 * ones(size(state));
end

function state = StateTransition(state, time)
 % rescale from [0.5, 2.5] to [-60, 60]
    minSet = 0.5 * ones(size(state));
    maxSet = 2.5 * ones(size(state));
    state = 120*(state - minSet) ./ (maxSet - minSet) - 60*ones(size(state));
    state = cappe_benchmark.StateTransition(state, time);
    % rescale from [-60, 60] to [0.5, 2.5]
    minSet2 = -60 * ones(size(state));
    maxSet2 = 60 * ones(size(state));
    state = (2.5 - 0.5) * (state - minSet2) ./ (maxSet2 - minSet2) + 0.5 * ones(size(state));
end
function y = GetDataC(x)
    y = sin(10*pi*x)./(2*x) + (x - 1).^4;
end

function y = GetData(x)
    y = sin(10*pi*x)./(2*x) + (x - 1).^4 + normrnd(0, 0.05, size(x));
end