classdef gpml_sode_wrapper < handle
    %GPML_SODE_WRAPPER Summary of this class goes here
    %   This is a convenience class to enable cleaner codes when doing
    %   sequential design of experiments work.  Not a necessary class though.
    
    %   This class is designed to contain useful functions for modifying
    %   a Gaussian Process model that was constructed using GPML
    %   Disclaimer, GPML is fairly flexible, and I have no idea what
    %   types of mean functions or covariance functions may break this
    %   wrapper.
    
    properties
        trainInput                                                          % Matrix: nSamples x nDimensions - contains the training data from the input space
        trainOutput                                                         % Matrix: nSamples x 1 - contains the training data from the output space
        testInput                                                           % Matrix: nTestSamples x nDimensions - contains points we might want to test the GP at
        testOutput                                                          % Matrix: nTestSamples x 1 - contains the output of the GP at testOutput points
        gpFuncHandle                                                        % Handle: function handle to evaluate the gaussian process, which must take three arguments (xx, yy, zz) See test1DSequentialDOE.m for example
        inputSampler                                                        % Handle: function that takes in one argument (number of samples) and returns random samples from the input space
        covfunc
        hyperparams
        K
        Cinv
        jitter
    end
    
    methods
        function obj = gpml_sode_wrapper(gpFuncHandle, trainInput,...
            trainOutput, testInput, inputSampler, covfunc, hyperparams)
            %GPML_SODE_WRAPPER Construct an instance of this class
            %   Detailed explanation goes here
            obj.gpFuncHandle = gpFuncHandle;
            obj.trainInput = trainInput;
            obj.trainOutput = trainOutput;
            obj.testInput = testInput;
            obj.testOutput = obj.predict(testInput);
            obj.inputSampler = inputSampler;
            
            obj.covfunc = covfunc;
            obj.hyperparams = hyperparams;
            [nTrain, ~] = size(trainInput);
            obj.jitter = exp(hyperparams.lik)*eye(nTrain);
            obj.K = obj.computeCovarianceMatrix() + obj.jitter;
            obj.Cinv = (obj.K)^-1;
        end
        
        function [meanPrediction, stdPrediction] = predict(obj, inputArray)
            %% predict
            % Uses the current training data to make a single prediction
            % from the gaussian process for a particular input array
            [meanPrediction, stdPrediction] = obj.gpFuncHandle(...
                obj.trainInput, obj.trainOutput, inputArray);
        end
        function [meanPrediction, stdPrediction] = predictUpdated(obj, inputArray,...
            newTrainIn, newTrainOut)
            %% predictUpdated
            % Uses the current training data AND new training data to make
            % predictions without modifying the base object.
            % Allows for rapidly testing out what new training point to add.
            %% NOTE: THIS WILL BE DEPRECATED AS A FASTER FUNCTION IS CORRECTLY CODED
            [meanPrediction, stdPrediction] = obj.gpFuncHandle(...
                [obj.trainInput; newTrainIn], [obj.trainOutput; newTrainOut],...
                inputArray);
        end
        
        function obj = updateFuncHandle(obj, newHandle)
            %% updateFuncHandle
            % b/c we do not store the parameters of the GP (hyper params,
            % mean functions, etc.,) we cannot do minimization in this 
            % class wrapper.  Thus, optimization of the hyperparameters
            % must be done externally, and then the function that this 
            % wrapper uses must be updated.
            obj.gpFuncHandle = newHandle;
        end
        
        function gpSamples = sampleGP(obj, nSamples)
            %% Generates samples of the output from the Gaussian process 
            % model using the original training data.
            inSamples = obj.inputSampler(nSamples);
            [tempMu, tempStd] = obj.predict(inSamples);
            gpSamples = normrnd(tempMu, tempStd);
        end
        
        function gpSamples = sampleUpdatedGP(obj, nSamples, ...
                newTrainIn, newTrainOut)
            %% Generates samples of the output from the Gaussian process 
            % model using updated training data.
            inSamples = obj.inputSampler(nSamples);
            [tempMu, tempStd] = obj.predictUpdated(inSamples, ...
                newTrainIn, newTrainOut);
            gpSamples = normrnd(tempMu, tempStd);
        end
        
        function pd = computePDF(obj, output, design, nDesign)
            % nDesign is the number of points in the design matrix
            [mus, stds] = obj.predict(design);                                  % Returns a mean and standard deviation for each input in the design matrix            
            pd = sum(normpdf(output, mus, stds))/nDesign;                       % Compute the PDF of each normal distribution evaluated at the output value
        end
        
        function pd = computeUpdatedPDF(obj, output, design, ...
                newTrainIn, newTrainOut, nDesign)
            % nDesign is the number of points in the design matrix
            [mus, stds] = obj.predictUpdated(...
                design, newTrainIn, newTrainOut);                           % Returns a mean and standard deviation for each input in the design matrix          
            pd = sum(normpdf(output, mus, stds)) / nDesign;                 % Compute the PDF of each normal distribution evaluated at the output value
        end
        function [stdDiff, stdDiffErr] = activeLearn(obj, alternatives, testPoints)
            [~, prevStd] = obj.predict(testPoints);
            [nAlternatives, ~] = size(alternatives);
            [nTest, ~] = size(testPoints);
            stdDiff = zeros([nAlternatives, nTest]);
            for iAlternative = 1:nAlternatives
                [~, tempStd] = obj.predictUpdated(testPoints, ...
                    alternatives(iAlternative, :), 1);
                stdDiff(iAlternative, :) = prevStd - tempStd;
            end
            stdDiffErr = std(stdDiff, 1, 2);
            stdDiff = mean(stdDiff, 2);
        end
        function K = computeCovarianceMatrix(obj)
            K = feval(obj.covfunc{:}, obj.hyperparams.cov, obj.trainInput) + obj.jitter;
        end
        function Kz = computeCovarianceVector(obj, input)
            Kz = feval(obj.covfunc{:}, obj.hyperparams.cov, obj.trainInput, input);
        end
        function Kx = computeSelfCovariance(obj, input)
            Kx = feval(obj.covfunc{:}, obj.hyperparams.cov, input, 'diag');
        end
        function KK = computeSoloCovarianceMatrix(obj, input1, input2)
            KK = feval(obj.covfunc{:}, obj.hyperparams.cov, input1, input2);
        end
        function [klDivArray, newOutputs] = KL_divergence_MonteCarlo(...
                obj, inputAlternatives, nSamples, nMC_samples, nDesign)
            % Estimates the KL_divergence between two Gaussian process models
            % using Monte Carlo Sampling
            %% Generate a set of possible updated Gaussian Process models
            % There are nAlternative possible inputs we could add (i.e., evaluate our high fidelity model),
            % and there are nSamples possible outputs that we could obtain if we evaluated our high fidelity
            % model for a particular input.  Thus, for our Monte Carlo simulation, we compute the KL divergence
            % for nAlternatives x nSamples possible updated Gaussian Process models.

            [nAlternatives, ~] = size(inputAlternatives);
            klDivArray = zeros(nAlternatives, nSamples);
            
            for iAlternative = 1:nAlternatives
                [mu, std] = obj.predict(inputAlternatives(iAlternative, :));               
                newOutputs = normrnd(mu, std, nSamples, 1);                 % Input used to create the new GP
                for iSample = 1:nSamples
                    fprintf("Computing KL Divergence for alternative %d : Output %d \n", iAlternative, iSample)
                    klDivArray(iAlternative, iSample) = obj.computeAverageKL(...
                        inputAlternatives(iAlternative,:),...
                        newOutputs(iSample), ...
                        nMC_samples, nDesign);
                end
            end
        end
        
        function avg_kl_div = computeAverageKL(obj, newTrainIn, ...
                newTrainOut, nOutSamples, nDesign)
            design = obj.inputSampler(nDesign);                             % Evaluates MC estimate of the PDF function
            
            gpSamples = obj.sampleUpdatedGP(nOutSamples, ...                % Points to evaluate the PDF function
                        newTrainIn, newTrainOut); 
            avg_kl_div = 0;
            for iOutput = 1:nOutSamples
                newPDF = obj.computeUpdatedPDF(gpSamples(iOutput),...
                    design, newTrainIn, newTrainOut, nDesign);
                oldPDF = obj.computePDF(gpSamples(iOutput), design, nDesign);
                avg_kl_div = avg_kl_div + log(newPDF) - log(oldPDF);
            end
            avg_kl_div = avg_kl_div / nOutSamples;
        end
    end
end
