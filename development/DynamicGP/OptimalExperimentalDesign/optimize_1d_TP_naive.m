clear all
clc
close all

%% User selected parameters
% parameters for controlling capability model
constantLoad = 4;
gpParamFile = 'dataFiles/base_line_gp_2.mat';
trainingInputTemp = linspace(0, 1, 5); % Initial Training Inputs for GP
trainingInputs = getStrains(trainingInputTemp, constantLoad*ones(size(trainingInputTemp)));
trainingOutputs = getCapability(trainingInputTemp); % Initial Training Outputs for GP
%testInputs = getStrains(normrnd(0.6, 0.05, 20), constantLoad*ones(20, 1)); % Normal Prior on damage
testInputs = getStrains(linspace(0, 1, 20), constantLoad*ones(20,1)); % Uniform Prior on damage
numTrials3 = 1000;
noiseLevel = 1e-2;
librarySize = 10;
% parameters for controlling knowledge gradient optimization
numTrials = 10;
numTrials2 = 100;
numExperiments = 50;
%array1 = linspace(0.13, 0.99, 9)';% training inputs
array2 = linspace(0, 1, 200)';% additional alternatives
damageArray = array2; % used to grab alternatives

%% Initialize Problem
% select a subset of the alternatives to train a Gaussian Process on
trainAlternatives = array2;
additionalAlternatives = array2;
alternatives = array2;
for idx = 1:length(alternatives(:,1))
    alternatives(idx, 1) = getStrains(alternatives(idx,1), constantLoad(1));
end
for idx = 1:length(trainAlternatives(:,1))
    trainAlternatives(idx, 1) = getStrains(trainAlternatives(idx,1), constantLoad(1));
end



%% Plot original GP
params = load(gpParamFile);
hyp = params.hyp;
meanfunc = params.meanfunc;
covfunc = params.covfunc;
likfunc = params.likfunc;

M = length(alternatives);
entropies = zeros(M, numTrials);
  hyp = minimize(hyp, @gp, -100, @infGaussLik,...
    meanfunc, covfunc, likfunc, trainingInputs, trainingOutputs);

hyp.lik(1) = -2.5;
testDam = linspace(0, 1, 200)';
test = getStrains(...
    testDam, constantLoad*ones(size(testDam)));
[randomMean, randomStdDev] = gp(...
                                            hyp, @infGaussLik, ...
                                            meanfunc, covfunc, ...
                                            likfunc, trainingInputs,...
                                            trainingOutputs, ...
                                            test...
                                            );

fig = figure('rend','painters','pos',[10 10 600 971], 'color', 'white');
subplot1 = subplot(2,1,1);

errorbar(test', randomMean', 2*randomStdDev')
xlim([0.18, 0.42])
ylabel('Capability', 'fontsize', 16)
xlabel('Strain', 'fontsize', 16)
hold on
%scatter(testInputs, zeros(length(testInputs), 1), 'red')
hold on


%% Offline Library Optimization Loop

frameCount = 1;

    % Build a Gaussian Process for expected entropy gain
trackLibrary = [];

for zz = 1:librarySize - length(trainingInputs)
    subplot1 = subplot(2, 1, 1);
    [maxi, index] = max(abs(randomStdDev));
    addIndex = 40;
    if zz == 1
        scatter(test(index+addIndex), randomMean(index+addIndex), 'red', 'filled')
    end
    if zz > 1
        cla(subplot2)
    end
    subplot2 = subplot(2, 1, 2);
    errorbar(test', randomMean', 2*randomStdDev')
    ylabel('Capability', 'fontsize', 16)
    xlabel('Strain', 'fontsize', 16)
    xlim([0.18, 0.42])
    hold on
    [newLibInput, newLibOutput] = updateGP(testDam(index+addIndex), ...
        constantLoad);
    
    trainingInputs = [trainingInputs; newLibInput];
    trainingOutputs = [trainingOutputs; newLibOutput];
    [randomMean, randomStdDev] = gp(...
                                            hyp, @infGaussLik, ...
                                            meanfunc, covfunc, ...
                                            likfunc, trainingInputs,...
                                            trainingOutputs, ...
                                            test...
                                            );
    drawnow
    frames(frameCount) = getframe(fig);
    frameCount = frameCount + 1;
    waitforbuttonpress
end


function [strains, capability] = updateGP(damage, load)
    strains = getStrains(damage, load);
    capability = getCapability(damage);
end