%% testDynamicGP_filteringCost.m
% Author: Brian Burrows
% This method is to test the cost, and accuracy of various surrogates
% using a fixed training set size

clear all
close all
clc
nRepeats = 100;
for iRepeat = 1:nRepeats

%% Problem Initialization
nUpdates = 5;
nClusters = 10;


lowerBound = -40;
upperBound = 40;
librarySize = 200;
nMonteCarlo = 100;
Tmax = 70;
libraryIn = linspace(lowerBound, upperBound, librarySize)';                  % Fix the
libraryOut = cappe_benchmark.GetDataC(libraryIn);

%% Filtering Initialization for Cappe Test Function
initialState = 0;
initialUncertainty = 10;
processNoise = 10;
measurementNoise = 1;
initialMeasurement = cappe_benchmark.GetData(initialState);

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
noisyOutputs = cappe_benchmark.GetData(libraryIn')';
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
trainIn =  libraryIn(1:25:end);                                 % only place the bounds of the problem in the training set
trainOut = libraryOut(1:25:end);                                % only place the outputs corresponding to the bounds of the training set

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


%% Initiazile Bootstrap Particle Filter
numParticles = 500;                                                        % Set the number of particles used at each iteration.  Increase for better accuracy but higher computational cost
pdf_g = @(dat,xx) normpdf(dat, 0.05 * xx^2, 1);                             % Set up the likelihood model

PF1.resampling = 1;                                                          % Set resampling = true, this will change dynamically later
PF1.likelihood = pdf_g; 
PF1.proposalSampler.random = @(xx, dat, Tim) normrnd(...                     % Create a random number generator to sample from the state evolution function
    cappe_benchmark.StateTransitionC(xx, Tim),...
    processNoise^0.5);
PF1.particles = zeros(length(initialState), numParticles);                   
PF1.weights = zeros(1, numParticles);
for idx = 1:numParticles
    data = cappe_benchmark.GetData(initialState);                           %
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
    PF2.likelihood = @(data, possibleStates) dynGP2.likelihood(data, possibleStates);
    PF3.likelihood = @(data, possibleStates) dynGP3.likelihood(data, possibleStates);
    

    % Evaluate the Bootstrap Particle Filter
    PF1.data = measurements(iTime - 1);
    PF1.oper = iTime - 1;
    tic
    PF1 = BootstrapFilter(PF1);
    pf1Time(iTime - 1) = toc;
    
    pf1Est(iTime - 1) = PF1.estimate;
    pf1Err(iTime - 1) = PF1.err^0.5;
    
    % Evaluate the Bootstrap Particle Filter
    PF2.data = measurements(iTime - 1);
    PF2.oper = iTime - 1;
    tic
    PF2 = BootstrapFilter(PF2);
    pf2Time(iTime - 1) = toc;
    
    pf2Est(iTime - 1) = PF2.estimate;
    pf2Err(iTime - 1) = PF2.err^0.5;
    
    % Evaluate the Bootstrap Particle Filter
    PF3.data = measurements(iTime - 1);
    PF3.oper = iTime - 1;
    tic
    PF3 = BootstrapFilter(PF3);
    pf3Time(iTime - 1) = toc;
    
    pf3Est(iTime - 1) = PF3.estimate;
    pf3Err(iTime - 1) = PF3.err^0.5;

    
    %% Evolve state and collect measurement
    if iTime < Tmax
        states(iTime) = cappe_benchmark.StateTransition(states(iTime -1), ...
            iTime);
        measurements(iTime) = cappe_benchmark.GetData(states(iTime));
    end
    %% Evolve the Gaussian Process Model to be better at the next time step
    tic
    dynGP2.inputSampler = @(numSampl) cappe_benchmark.StateTransition(...
        PF2.particles(:, 1:numSampl), iTime)';
    dynGP2.updateGaussianProcess(numParticles, nUpdates);
    updateTime2(iTime) = toc;
    
    tic
    dynGP3.inputSampler = @(numSampl) cappe_benchmark.StateTransition(...
        PF3.particles(:, 1:numSampl), iTime)';
    dynGP3.updateGaussianProcess_kMeans(numParticles, nUpdates, nClusters);
    updateTime3(iTime) = toc;
    
    fprintf("Updating Gaussian Process took %f seconds \n", toc)
    toc
end


e1(iRepeat) = mean((states - pf1Est).^2)^0.5;
e2(iRepeat) = mean((states - pf2Est).^2)^0.5;
e3(iRepeat) = mean((states - pf3Est).^2)^0.5;

end

figure(1)
subplot(3,1,1)
niceErrorbar([1:Tmax]', zeros([Tmax,1]), pf1Err')
hold on
scatter([1:Tmax], states - pf1Est, 10, 'filled', 'k');
plot([1:Tmax], e1*ones([1, Tmax]))

subplot(3,1,2)
niceErrorbar([1:Tmax]', zeros([Tmax,1]), pf2Err')
hold on
scatter([1:Tmax], states - pf2Est, 10, 'filled', 'k');
plot([1:Tmax], e2*ones([1, Tmax]))

subplot(3,1,3)
niceErrorbar([1:Tmax]', zeros([Tmax,1]), pf3Err')
hold on
scatter([1:Tmax], states - pf3Est, 10, 'filled', 'k');
plot([1:Tmax], e3*ones([1, Tmax]))
