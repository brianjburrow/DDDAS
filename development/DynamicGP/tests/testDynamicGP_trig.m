%% testDynamicGP.m
% Author: Brian Burrows
% This class is to demonstrate the functionality of the Dynamic GP class
% and to make sure that new functions in the class is working appropriately

clear all
close all
clc

%% Problem Initialization
nUpdates = 5;
nClusters = 10;
    
v = VideoWriter('dynamicGP_filtering.avi');
open(v);

lowerBound = -40;
upperBound = 40;
librarySize = 200;
maxTrainingSize = 15;
nMonteCarlo = 100;
Tmax = 70;
libraryIn = linspace(lowerBound, upperBound, librarySize)';                  % Fix the
libraryOut = cappe_benchmark_trig.GetDataC(libraryIn);

%% Filtering Initialization for Cappe Test Function
initialState = 0;
initialUncertainty = 10;
processNoise = 10;
measurementNoise = 1;

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
noisyOutputs = cappe_benchmark_trig.GetData(libraryIn')';
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
trainIn =  libraryIn(1:20:end);                                 % only place the bounds of the problem in the training set
trainOut = libraryOut(1:20:end);                                % only place the outputs corresponding to the bounds of the training set

dynGP = dynamicGP(gpHandle, trainIn,...
            trainOut, inputSampler,...
            libraryIn, libraryOut);

%% Train the initial GP
disp(maxTrainingSize - length(trainIn))
dynGP = dynGP.addNewPoints(nMonteCarlo, maxTrainingSize - length(trainIn));
[muVec, stdVec] = dynGP.predict(libraryIn);

figure(1)
subplot(3, 2, 1)
niceErrorbar(libraryIn, muVec, 3*stdVec)
hold on
ylim([-10, 100])
xlim([lowerBound, upperBound])
title("Current GP")
plot(libraryIn, libraryOut, 'k', 'LineStyle', '-.', 'LineWidth', 2.0)
scatter(dynGP.trainInput, dynGP.trainOutput, 20, 'b', 'filled')
hold off

subplot(3,2,3)
title("Filtering Distribution")
hold on
area(libraryIn, normpdf(libraryIn, initialState, initialUncertainty^0.5), 'FaceColor', 'blue');
xlim([lowerBound, upperBound])
hold off

subplot(3,2,5)
xlabel('True State')
hold on
ylabel('True Data')
plot(libraryIn, libraryOut, 'k', 'LineStyle', '-', 'LineWidth', 0.5)
xlim([lowerBound, upperBound])
ylim([-10, 100])
initialMeasurement = cappe_benchmark_trig.GetData(initialState);
scatter(initialState, initialMeasurement, 20, 'r', 'filled')
hold off

%% Initiazile Bootstrap Particle Filter
numParticles = 2000;                                                        % Set the number of particles used at each iteration.  Increase for better accuracy but higher computational cost
pdf_g = @(dat,xx) normpdf(dat, 0.05 * xx^2, 1);                             % Set up the likelihood model

PF.resampling = 1;                                                          % Set resampling = true, this will change dynamically later
PF.likelihood = pdf_g; 
PF.proposalSampler.random = @(xx, dat, Tim) normrnd(...                     % Create a random number generator to sample from the state evolution function
    cappe_benchmark_trig.StateTransitionC(xx, Tim),...
    processNoise^0.5);
PF.particles = zeros(length(initialState), numParticles);                   
PF.weights = zeros(1, numParticles);
for idx = 1:numParticles
    data = cappe_benchmark_trig.GetData(initialState);                           %
    PF.particles(:,idx) = PF.proposalSampler.random(initialState, ...       % generate an initial particle
        data, 0);
    PF.weights(idx) = PF.likelihood(data, initialState);                    % Compute initial weights
end
PF.weights = PF.weights/sum(PF.weights);                                    % Normalize the weights

states = initialState;
errs = initialUncertainty^0.5;
measurements = initialMeasurement;
for iTime = 2:Tmax
    fprintf("Processing Time Step %d \n", iTime)
    %% Perform Particle filtering
    PF.likelihood = @(data, possibleStates) dynGP.likelihood(data, possibleStates')';
    
    % Evaluate the Bootstrap Particle Filter
    PF.data = measurements(iTime - 1);
    PF.oper = iTime;
    tic
    PF = BootstrapFilter(PF);
    pfTime(iTime) = toc;
    
    pfEst(iTime) = PF.estimate;
    pfErr(iTime) = PF.err^0.5;
    
    %% After prediction, update the pdf of the state
    if iTime > 2
        % resample particles for plotting
        tempPart = zeros(size(PF.particles));
        count = 1;
        r = mnrnd(numParticles, PF.weights);                                % Multinomial resampling strategy (Not recommended in general)
                                                                            % May be useful to replace this with some other strategy
        for idx = 1:numParticles
            for dmx = 1:r(idx)
                tempPart(:, count) = PF.particles(:,idx);
                count = count + 1;
            end
        end

        subplot(3, 2, 1)
        niceErrorbar(libraryIn, tempMU, 3*tempStd)
        hold on
        ylim([-10, 100])
        xlim([lowerBound, upperBound])
        title("Current GP")
        plot(libraryIn, libraryOut, 'k', 'LineStyle', '-.', 'LineWidth', 2.0)
        scatter(dynGP.trainInput, dynGP.trainOutput, 20, 'b', 'filled')
        hold off
        
        subplot(3, 2, 3);
        histogram(tempPart)
        hold on
        xlim([lowerBound, upperBound])
        title("Filtering Distribution")
        hold off
        
        
        subplot(3,2,5)
        plot(libraryIn, libraryOut, 'k', 'LineStyle', '-', 'LineWidth', 0.5)
        hold on
        xlabel("True State")
        hold on
        ylabel("True Data")
        xlim([lowerBound, upperBound])
        ylim([-10, 100])
        initialMeasurement = cappe_benchmark_trig.GetData(initialState);
        scatter(states(iTime - 1), measurements(iTime - 1), 20, 'r', 'filled')
        hold off
        
        frame = getframe(gcf);
        for ii = 1:10
            writeVideo(v,frame);
        end
    end
    
    %% Evolve state and collect measurement
    states(iTime) = cappe_benchmark_trig.StateTransition(states(iTime -1), ...
        iTime);
    measurements(iTime) = cappe_benchmark_trig.GetData(states(iTime));
    
    %% Evolve the Gaussian Process Model to be better at the next time step
    tic
    dynGP.inputSampler = @(numSampl) cappe_benchmark_trig.StateTransition(...
        PF.particles(:, 1:numSampl), iTime)';
    %dynGP.updateGaussianProcess(numParticles, nUpdates);
    dynGP.updateGaussianProcess_kMeans(numParticles, nUpdates, nClusters);
    fprintf("Updating Gaussian Process took %f seconds \n", toc)
    toc
    
    %% Plotting
    subplot(3, 2, 4);
    histogram(dynGP.inputSampler(numParticles))
    hold on
    xlim([lowerBound, upperBound])
    title("Forecast Distribution")
    hold off
    subplot(3, 2, 2)
    [tempMU, tempStd] = dynGP.predict(libraryIn);
    niceErrorbar(libraryIn, tempMU, 3*tempStd)
    hold on
    plot(libraryIn, libraryOut, 'k', 'LineStyle', '-.', 'LineWidth', 2.0)
    scatter(dynGP.trainInput, dynGP.trainOutput, 20, 'b', 'filled')
    xlim([lowerBound, upperBound])
    ylim([-20, 100])
    title("Forecasted GP")
    hold off
    
    subplot(3, 2, 6)
    tSet = [1:iTime]';
    eSet = [pfEst - states(1:iTime)]';
    pfSet = pfErr';
    niceErrorbar(tSet, zeros([iTime, 1]), 3*pfSet)
    hold on
    scatter(tSet, eSet, 15, 'k', 'filled');
    xlim([1, Tmax])
    ylim([-40, 40])
    title("Filtering Error")
    hold off
    drawnow
    frame = getframe(gcf);
    for ii = 1:10
        writeVideo(v,frame);
    end
    
end
close(v)