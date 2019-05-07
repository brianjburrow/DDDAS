classdef ActiveLearning < handle
    %ActiveLearning Summary of this class goes here
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
        
        function [meanPrediction, stdPrediction] = predictReduced(obj, inputArray,...
            index)
            %% predictUpdated
            % Uses the current training data AND new training data to make
            % predictions without modifying the base object.
            % Allows for rapidly testing out what new training point to add.
            %% NOTE: THIS WILL BE DEPRECATED AS A FASTER FUNCTION IS CORRECTLY CODED
            
            [nTrain, ~] = size(obj.trainInput);
            if index == 1
                tempTrainInput = obj.trainInput(2:end, :);
                tempTrainOutput = obj.trainOutput(2:end, :);
            elseif index == nTrain
                tempTrainInput = obj.trainInput(1:nTrain - 1, :);
                tempTrainOutput = obj.trainOutput(1:nTrain - 1, :);
            else
                tempTrainInput = [obj.trainInput(1:index - 1, :);
                    obj.trainInput(index + 1:end,:, :)];
                
                tempTrainOutput = [obj.trainOutput(1:index - 1, :);
                    obj.trainOutput(index + 1:end,:, :)];
            end

            [meanPrediction, stdPrediction] = obj.gpFuncHandle(...
                tempTrainInput, tempTrainOutput,...
                inputArray);
        end
        
        function obj = updateTrainingData(obj, newTrainIn, newTrainOut)
            obj.trainInput = [obj.trainInput; newTrainIn];
            obj.trainOutput = [obj.trainOutput; newTrainOut];
        end
        
        function obj = removeTrainingData(obj, index)
            obj.trainInput(index, :) = [];
            obj.trainOutput(index, :) = [];
        end
        
        function obj = updateFuncHandle(obj, newFuncHandle)
            %% updateFuncHandle
            % b/c we do not store the parameters of the GP (hyper params,
            % mean functions, etc.,) we cannot do minimization in this 
            % class wrapper.  Thus, optimization of the hyperparameters
            % must be done externally, and then the function that this 
            % wrapper uses must be updated.
            obj.gpFuncHandle = newFuncHandle;
        end
      
        function [selection, stdDiff, stdDiffErr] = activeLearn(obj, alternatives, nTest)
            testPoints = obj.inputSampler(nTest);
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
            [~, selection] = max(stdDiff);
        end
        
        function [selection, stdDiff, stdDiffErr] = activeUnlearn(obj, nTest)
            alternatives = obj.trainInput;
            testPoints = obj.inputSampler(nTest);
            [~, prevStd] = obj.predict(testPoints);
            [nAlternatives, ~] = size(alternatives);
            [nTest, ~] = size(testPoints);
            stdDiff = zeros([nAlternatives, nTest]);
            for iAlternative = 1:nAlternatives
                [~, tempStd] = obj.predictReduced(testPoints, ...
                    iAlternative);
                stdDiff(iAlternative, :) = tempStd - prevStd;
            end
            stdDiffErr = std(stdDiff, 1, 2);
            stdDiff = mean(stdDiff, 2);
            [~, selection] = min(stdDiff);
        end
        %% Not Super useful add ons
       function gpSamples = sampleGP(obj, nSamples)
            %% Generates samples of the output from the Gaussian process 
            % model using the original training data.
            inSamples = obj.inputSampler(nSamples);
            [tempMu, tempStd] = obj.predict(inSamples);
            gpSamples = normrnd(tempMu, tempStd);
       end
        
       function obj = sortTrainingData(obj)
           [~, nDim] = size(obj.trainInput);
           if nDim ~= 1
               error('sortTrainingData not implemented for input spaces with more than one dimension')
           end
           [obj.trainInput, obj.trainOutput] = sortTwoArrays(obj.trainInput, obj.trainOutput);
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
    end
end

