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
testInputs = getStrains(normrnd(0.6, 0.05, 20), constantLoad*ones(20, 1)); % Normal Prior on damage
%testInputs = getStrains(linspace(0, 1, 20), constantLoad*ones(20,1)); % Uniform Prior on damage

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
test = linspace(0.19, 0.41, 200)';
[randomMean, randomStdDev] = gp(...
                                            hyp, @infGaussLik, ...
                                            meanfunc, covfunc, ...
                                            likfunc, trainingInputs,...
                                            trainingOutputs, ...
                                            test...
                                            );

fig = figure('rend','painters','pos',[10 10 600 971], 'color', 'white');
subplot1 = subplot(3,1,1);

errorbar(test', randomMean', 2*randomStdDev')
xlim([0.18, 0.42])
ylabel('Capability')
hold on
%scatter(testInputs, zeros(length(testInputs), 1), 'red')
hold on

subplot3 = subplot(3, 1, 3);
errorbar(test', randomMean', 2*randomStdDev')
ylabel('Capability')
xlabel('Strain')
xlim([0.18, 0.42])
%% Offline Library Optimization Loop

frameCount = 1;

    % Build a Gaussian Process for expected entropy gain
trackLibrary = [];
trackLibraryOutput = [];
for zz = 1:librarySize - length(trainingInputs)
    %% Estimate Sampling Statistics
    [means, samplingVars, stdError]  = estimateSamplingStatistics(...
    trainAlternatives, numTrials, numTrials2, numTrials3, gpParamFile, ...
    trainingInputs, trainingOutputs, testInputs, noiseLevel);
    means = means';
   
    hyp4 = struct('mean', [], ...
        'cov', log([0.01, 2*max(means).^0.5]),...
        'lik', log(max(abs(stdError))));
    meanfunc2 = {};
    covfunc2 = {@covSEiso};
    likfunc2 = @likGauss;

    hyp2 = minimize(hyp4, @gp, -100, @infGaussLik,...
        meanfunc2, covfunc2, likfunc2, trainAlternatives, means);
    hyp2.cov(1) = log(0.01);
    hyp2.lik(1) = log(3*max(abs(stdError)));
    %disp(log(max(stdError)));

    [means2, ~] = gp(hyp2, @infGaussLik, meanfunc2, covfunc2, likfunc2, ...
        trainAlternatives, means, alternatives);

    % get full covariance matrix
    [negloglik, derivnegloglik, poststr] = gp(...
        hyp2, @infGaussLik, meanfunc2, covfunc2, likfunc2, ...
        trainAlternatives, means);
    kTT = feval(covfunc2{:}, hyp2.cov, alternatives);
    kTI = feval(covfunc2{:}, hyp2.cov, alternatives, trainAlternatives);
    kII = feval(covfunc2{:}, hyp2.cov, trainAlternatives);

    ssig = exp(2*hyp2.lik)*eye(max(size(trainAlternatives)));
    covarianceMatrix = kTT - kTI/(kII + ssig)*kTI';
    
    
  


    % Inner Loop: Knowledge Gradient Optimization
    %errorbar(alternatives, means2, samplingVars);

    track_experiments = [];
    track_experiment_output = [];
    subplot2 = subplot(3,1,2);
    
    for kk = 1:numExperiments
        [x_star, nu, a] = CorrelatedKG(means2, covarianceMatrix,...
            samplingVars);
        % new experiment
        [new_data, ~, ~] = estimateSamplingStatistics(alternatives(x_star),...
            numTrials, numTrials2, numTrials3, gpParamFile, trainingInputs, ...
            trainingOutputs, testInputs, noiseLevel);
        means2 = updateMean(means2, covarianceMatrix, new_data, x_star, ...
            samplingVars);
        covarianceMatrix = updateVar(covarianceMatrix, x_star, samplingVars);
       
        cla(subplot2)

        disp('plotting')
        hold on 
        
        errs = 3*diag(covarianceMatrix).^0.5;
        [tempIn, sortIndex] = sort(alternatives);
        errorbar(alternatives(sortIndex), means2(sortIndex), ...
            errs(sortIndex));
        plot(alternatives(sortIndex), means2(sortIndex), 'black')
        ylabel('K-L Divergence')
        xlim([0, 1])
        hold on
 
        if kk == 1
            track_experiments(kk) = x_star;
            track_experiment_output(kk) = new_data;
            scatter(alternatives(track_experiments), new_data, 'red', 'filled')
        else
            scatter(alternatives(track_experiments), track_experiment_output, 'black','filled')
            track_experiments(kk) = x_star;
            track_experiment_output(kk) = new_data;
            scatter(alternatives(x_star), new_data, 'red', 'filled')
        end
        xlim([0.18, 0.42])
        frames(frameCount) = getframe(fig);
        frameCount = frameCount + 1;
    end
    [maxKLdiv, indexToAdd] = max(means2);
    disp(alternatives(indexToAdd,:))
    [newLibInput, newLibOutput] = updateGP(damageArray(indexToAdd), ...
        constantLoad);
    
    scatter(alternatives(indexToAdd), maxKLdiv, 'green', 'filled')
    trainingInputs = [trainingInputs; newLibInput];
    trainingOutputs = [trainingOutputs; newLibOutput];
    
    subplot3 = subplot(3,1,3);
    cla(subplot3)
 
    [randomMean, randomStdDev] = gp(...
                                            hyp, @infGaussLik, ...
                                            meanfunc, covfunc, ...
                                            likfunc, trainingInputs,...
                                            trainingOutputs, ...
                                            test...
                                            );

                                       
    errorbar(test', randomMean', 2*randomStdDev')
    xlim([0.18, 0.42])
    hold on
    if zz == 1
        trackLibrary(zz) = newLibInput;
        trackLibraryOutput(zz) = newLibOutput;
        scatter(trackLibrary, trackLibraryOutput, 'green', 'filled')
    else
        scatter(trackLibrary, trackLibraryOutput, 'black', 'filled')
        trackLibrary(zz) = newLibInput;
        trackLibraryOutput(zz) = newLibOutput;
        scatter(newLibInput, newLibOutput, 'green', 'filled')
    end
    xlim([0.18, 0.42])
    ylabel('Capability')
    xlabel('Strain')
    for abc = 1:10
        frames(frameCount) = getframe(fig);
        frameCount = frameCount + 1;
    end

end


function [strains, capability] = updateGP(damage, load)
    strains = getStrains(damage, load);
    capability = getCapability(damage);
end