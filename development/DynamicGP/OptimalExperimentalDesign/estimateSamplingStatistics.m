function [means, samplingVars, stdError] = estimateSamplingStatistics(...
    alternatives, ...
    numTrials, ...
    numTrials2, ...
    numTrials3, ...
    gaussianProcessParamFile, ...
    trainingInputs, ...
    trainingOutputs, ...
    testInputs, ...
    noiseLevel ...
)

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
            
%             [whuh, whahh] = gp(...
%                 hyp, @infGaussLik, meanfunc, covfunc, likfunc, ...
%                 newTrainInput, newTrainOutput, test);
%             clf
%             ylim([0, 500])
%             errorbar(test, whuh, whahh)
%             ylim([0, 500])
%             drawnow
            [mu_new, sig_new] = convert_GMM_toGaussian(sup, ...
                                                     mu_array, ....
                                                     sig_array);

            MC_integration_samples = normrnd(mu_old, sig_old, numTrials);
            for zzz = 1:numTrials2
                kl_div = kl_div + log(...
                    normpdf(MC_integration_samples(zzz), mu_old, sig_old)/...
                    normpdf(MC_integration_samples(zzz), mu_new, sig_new));
            end
            kl_div = kl_div / numTrials2;
            entropies(idx, dmx) = kl_div;
        end
    end

means = mean(abs(entropies), 2)';
disp(size(entropies))
samplingVars = var(entropies')';
%stdError = (samplingVars).^0.5;
stdError = (samplingVars.^0.5) / (numTrials ^ 0.5);
end