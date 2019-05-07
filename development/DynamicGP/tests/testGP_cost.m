%% testGP_cost.m
% Author: Brian Burrows
% Test the cost of evaluating a Gaussian process model
% given changes to the training set size, the number of samples evaluated
% at once

clear all
close all
clc


%% Problem Initialization
nUpdates = 3;
nClusters = 20;


lowerBound = -40;
upperBound = 40;
librarySize = 100;
nMonteCarlo = 100;
Tmax = 70;
libraryIn = linspace(lowerBound, upperBound, librarySize)';                  % Fix the
libraryOut = cappe_benchmark.GetDataC(libraryIn);


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
trainIn =  libraryIn(1:50:end);                                 % only place the bounds of the problem in the training set
trainOut = libraryOut(1:50:end);                                % only place the outputs corresponding to the bounds of the training set

[maxTrainingSize, ~] = size(trainIn);
predictSizes = [100, 250, 500, 750, 1000];
trainSizes = 1:4:1000;
nPredicts = length(predictSizes);
nTrains = length(trainSizes);
times = zeros(nPredicts, nTrains);
err = zeros(nPredicts, nTrains);
numTrain = zeros(nPredicts, nTrains);
for iPredictSize = 1:length(predictSizes)
    disp(iPredictSize)
    predictSize = predictSizes(iPredictSize);
    for iTrainSize = 1:length(trainSizes)
        temp = zeros(1, 100);
        trainIn =  rand(trainSizes(iTrainSize),1);                                 % only place the bounds of the problem in the training set
        trainOut = rand(trainSizes(iTrainSize),1);                                % only place the outputs corresponding to the bounds of the training set
        dynGP1 = dynamicGP(gpHandle, trainIn,...
                    trainOut, inputSampler,...
                    libraryIn, libraryOut);
        for iRep = 1:100
            tic;
            dynGP1.predict(rand(predictSize, 1));
            temp(iRep) = toc;
        end
        times(iPredictSize, iTrainSize) = mean(temp);
        err(iPredictSize, iTrainSize) = std(temp);
        clear temp
        numTrain(iPredictSize, iTrainSize) = length(trainIn);
    end
end
plot(numTrain(1,:), times(1,:), 'k', 'LineWidth', 2.0)
hold on
plot(numTrain(2,:), times(2,:), 'r', 'LineWidth', 2.0)
plot(numTrain(3,:), times(3,:), 'g','LineWidth', 2.0)
plot(numTrain(4,:), times(4,:), 'b','LineWidth', 2.0)
plot(numTrain(5,:), times(5,:), 'm','LineWidth', 2.0)
legend('numPredict: 100', 'numPredict: 250', 'numPredict: 500', 'numPredict: 750', 'numPredict: 1000')
xlabel('Number of Training Points')
ylabel('Computational Cost')