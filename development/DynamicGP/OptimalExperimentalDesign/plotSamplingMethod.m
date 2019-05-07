clear all
clc
close all

%% User selected parameters
% parameters for controlling capability model
constantLoad = 4;
gpParamFile = 'dataFiles/base_line_gp_2.mat';
trainingInputTemp = linspace(0, 1, 10); % Initial Training Inputs for GP
trainingInputs = getStrains(trainingInputTemp, constantLoad*ones(size(trainingInputTemp)));
trainingOutputs = getCapability(trainingInputTemp); % Initial Training Outputs for GP
testInputs = getStrains(normrnd(0.6, 0.05, 20), constantLoad*ones(20, 1)); % Normal Prior on damage
%testInputs = getStrains(linspace(0, 1, 20), constantLoad*ones(20,1)); % Uniform Prior on damage
numTrials3 = 1000;
noiseLevel = 1e-2;
librarySize = 10;
% parameters for controlling knowledge gradient optimization
numTrials = 20;
numTrials2 = 10;
numExperiments = 50;
%array1 = linspace(0.13, 0.99, 9)';% training inputs
array2 = linspace(0, 1, 20)';% additional alternatives
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

%% Offline Library Optimization Loop

frameCount = 1;

[means, samplingVars, stdError, frameSet]  = plotSampling(...
    trainAlternatives, numTrials, numTrials2, numTrials3, gpParamFile, ...
    trainingInputs, trainingOutputs, testInputs, noiseLevel, array2);
    means = means';
disp(trainAlternatives)
aoeueu
function [means, samplingVars, stdError, frames] = plotSampling(...
    alternatives, ...
    numTrials, ...
    numTrials2, ...
    numTrials3, ...
    gaussianProcessParamFile, ...
    trainingInputs, ...
    trainingOutputs, ...
    testInputs, ...
    noiseLevel, ...
    damageAlternatives ...
)
    fig = figure('rend','painters','pos',[10 10 971 600], 'color', 'white');
    sp1 = subplot(1, 2, 1);
    % This function adds an alternative to the gaussian process
    % training data, and then evaluates the expected performance
    % of the new gaussian process on the test samples.
    % Inputs:
    %        alternatives: M x d list of inputs, M altenatives, d
    %                       dimensions
    %        numTrials   : scalar.  Number of monte carlo evaluations
    %                               used to evaluate the expectation
    %        gaussianProcess: GPML gaussian process object
    %        trainingInputs : N x d array.  N training points, d dimensions
    %        trainingOutputs: N x 1 array.
    %        testInputs: K x d array.  K test points, d dimensions
    %        gaussianProcessParamFile = filename where gp params are stored
    %                                   must have obj.hyp, obj.meanfunc,
    %                                   obj.likfunc.  Assumes @infGaussLik

    params = load(gaussianProcessParamFile);
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
 
    errorbar(test, randomMean, randomStdDev);
    drawnow
    
    % test original capability model on the test set
    [mu_array, sig_array] = gp(...
        hyp, @infGaussLik, meanfunc, covfunc, likfunc, ...
        trainingInputs, trainingOutputs, testInputs);
    % Generate GMM on resulting estimates, and convert back to 
    % Gaussian distribution
    sup = [];
    [mu_old,sig_old] = convert_GMM_toGaussian(sup, ...
                                                mu_array,...
                                                sig_array);
 
    prior_noisy_measurements = zeros(length(testInputs(:,1))*numTrials3, ...
        length(testInputs(1,:)));
    counter = 1;
    for idx = 1:length(testInputs(:,1))
        for dmx = 1:numTrials3
            prior_noisy_measurements(counter,:) = normrnd(...
                testInputs(idx,:), noiseLevel, 1); % generate a noisy strain
            counter = counter + 1;
        end
    end
   frameCount = 1;
    for idx = 1:M
        % get output of GP at alternative idx
        [randomMean, randomStdDev] = gp(...
                                                hyp, @infGaussLik, ...
                                                meanfunc, covfunc, ...
                                                likfunc, trainingInputs,...
                                                trainingOutputs, ...
                                                alternatives(idx,:)...
                                                );
        % sample a series of noisy samples from the GP output at alt idx
        newTrainInput = [trainingInputs; alternatives(idx,:)];
        noisy_output = normrnd(randomMean, randomStdDev, [numTrials, 1]);
        
        for dmx = 1:numTrials
            kl_div = 0;
            % append noisy output to training set
            newTrainOutput = [trainingOutputs; noisy_output(dmx)];
           
            % test updated GP on the test set                                   
            clear mu_array sig_array
            [mu_array, sig_array] = gp(...
                hyp, @infGaussLik, meanfunc, covfunc, likfunc, ...
                newTrainInput, newTrainOutput, testInputs);
            
            [whuh, whahh] = gp(...
                hyp, @infGaussLik, meanfunc, covfunc, likfunc, ...
                newTrainInput, newTrainOutput, test);
            cla(sp1)
            sp1 = subplot(1,2,1);
            ylim([0, 500])
            errorbar(test, whuh, whahh)
            ylim([0, 500])
            xlabel('Strain', 'fontsize', 18)
            ylabel('Capability', 'fontsize', 18)
            drawnow
            a = get(gcf,'Position');
            [mu_new, sig_new] = convert_GMM_toGaussian(sup, ...
                                                     mu_array, ....
                                                     sig_array);

            MC_integration_samples = normrnd(mu_old, sig_old, numTrials);
            for zzz = 1:numTrials2
                kl_div = kl_div + log(...
                    normpdf(MC_integration_samples(zzz), mu_old, sig_old)/...
                    normpdf(MC_integration_samples(zzz), mu_new, sig_new));
            end
            
            kl_div = kl_div/numTrials2;
            entropies(idx, dmx) = kl_div;
            if dmx == 1
                sp2 = subplot(1,2,2);
                scatter(damageAlternatives(idx), kl_div, 'red', 'filled');
                xlim([0, 1])
                hold on
                drawnow
            else
                sp2 = subplot(1,2,2);
                for zzz = idx:idx
                    scatter(damageAlternatives(zzz), ...
                        mean(entropies(zzz, 1:dmx)), 'red', 'filled'...
                        );
                    hold on
                    ylim([-0.0001,abs(max(max(entropies)))])
                    xlim([0, 1])
                    xlabel('Damage Alternative', 'fontsize', 18)
                    ylabel('Kl-Divergence', 'fontsize', 18)
                    drawnow

                    frames(frameCount) = getframe(fig);
                    frameCount = frameCount + 1;

                end
            end
        end
        
    end

means = mean(abs(entropies), 2)';
disp(size(entropies))
samplingVars = var(entropies')';
%stdError = (samplingVars).^0.5;
stdError = (samplingVars.^0.5) / (numTrials ^ 0.5);
end