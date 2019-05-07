%% testDynamicGP.m
% Author: Brian Burrows
% This class is to demonstrate the functionality of the Dynamic GP class
% and to make sure that new functions in the class is working appropriately

clear all
close all
clc

%% Problem Initialization
lowerBound = 0;
upperBound = 1;
librarySize = 200;
maxTrainingSize = 10;
nMonteCarlo = 100;

libraryIn = linspace(lowerBound, upperBound, librarySize)';                  % Fix the
libraryOut = testFuncDynGP(libraryIn);

%% Filtering Initialization
initialState = 0;
initialUncertainty = 10;

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
hyperParam = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                 covfunc, likfunc, libraryIn, libraryOut + normrnd(0, 0.1, [librarySize, 1])); 

%% Set up instance of the dynamicGP class
% set up GPML function handle
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, ...
    meanfunc, covfunc,...
    likfunc, xx,...
    yy, zz);
% Set up initial pdf over our input space
inputSampler = @(numSamp) normpdf(0.1, 0.1);
trainIn =  libraryIn(1:20:end);                                 % only place the bounds of the problem in the training set
trainOut = libraryOut(1:20:end);                                % only place the outputs corresponding to the bounds of the training set

dynGP = dynamicGP(gpHandle, trainIn,...
            trainOut, inputSampler,...
            libraryIn, libraryOut);

%% Train the initial GP
dynGP = dynGP.addNewPoints(nMonteCarlo, maxTrainingSize - length(trainIn));
[muVec, stdVec] = dynGP.predict(libraryIn);
figure(1)
subplot(2, 1, 1)
niceErrorbar(libraryIn, muVec, 3*stdVec)
hold on
title("Original GP")
plot(libraryIn, libraryOut, 'k', 'LineStyle', '-.', 'LineWidth', 2.0)
scatter(dynGP.trainInput, dynGP.trainOutput, 20, 'b', 'filled')
subplot(2,1,2)
plot(libraryIn, normpdf(libraryIn, 0.1, 0.1), 'b', 'LineStyle', '-.', 'LineWidth', 2.0);


function y = testFuncDynGP(x)
    y = sin(x*20);
end
