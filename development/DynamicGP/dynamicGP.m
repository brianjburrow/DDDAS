classdef dynamicGP < handle
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
        gpFuncHandle                                                        % Handle: function handle to evaluate the gaussian process, which must take three arguments (xx, yy, zz) See test1DSequentialDOE.m for example
        inputSampler                                                        % Handle: function that takes in one argument (number of samples) and returns random samples from the input space
        covfunc
        hyperparams
        K
        Cinv
        jitter
        libraryInputs
        libraryOutputs
    end
    methods        
        function obj = dynamicGP(gpFuncHandle, trainInput,...
            trainOutput, inputSampler, libraryInputs, libraryOutputs)
            % dynamicGP Constructor: Construct an instance of this class
            %   Detailed explanation goes here
            obj.libraryInputs = libraryInputs;
            obj.libraryOutputs = libraryOutputs;
            obj.gpFuncHandle = gpFuncHandle;
            obj.trainInput = trainInput;
            obj.trainOutput = trainOutput;
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
        
        function p = likelihood(obj, data, particles)
            [muPred, stdPred] = obj.predict(particles);
            p = normpdf(data, muPred, stdPred);
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
        function obj = updateGaussianProcess(obj, nTest, nChanges)
            alternatives = obj.libraryInputs;
            outputs = obj.libraryOutputs;
            %% Remove nChanges elements from the training data
            [selections, ~, ~] = obj.activeUnlearn(nTest, nChanges);
            obj.removeTrainingData(selections);
            for iChange = 1:nChanges
                [selections, ~, ~] = obj.activeLearn(alternatives, nTest);
                obj.updateTrainingData(alternatives(selections, :), outputs(selections, :));
                alternatives(selections, :) = [];
                outputs(selections, :) = [];
            end
        end
        
        function obj = updateGaussianProcess_kMeans(obj, nTest, nChanges, nClusters)
            particles = obj.inputSampler(nTest);
            %% Partition the number of particles into clusters
            tic
            [~, centroids] = kmeans(particles, nClusters, 'dist', 'sqeuclidean');
            
            [indx, ~] = knnsearch(obj.libraryInputs, centroids, 'K', 1);
            fprintf('reduction time %f', toc)
             indx = sort(indx);
             alternatives = obj.libraryInputs(indx,:);
             outputs = obj.libraryOutputs(indx,:);
              %% Remove nChanges elements from the training data
            [selections, ~, ~] = obj.activeUnlearn(nTest, nChanges);
            obj.removeTrainingData(selections);

            testPoints = obj.inputSampler(nTest);         
            [~, prevStd] = obj.predict(testPoints);
            [nAlternatives, ~] = size(alternatives);
            [nTest, ~] = size(testPoints);
            stdDiff = zeros([nAlternatives, nTest]);
            tempPrediction = @(alt) obj.predictUpdated(testPoints, alt, 1);
            for iChange = 1:nChanges
                for iAlternative = 1:nAlternatives
                    [~, tempStd] = tempPrediction(alternatives(iAlternative, :));
                    stdDiff(iAlternative, :) = prevStd - tempStd;
                end
                stdDiff2 = mean(stdDiff, 2); 
                [~, selection] = max(stdDiff2);
                
                obj.trainInput = [obj.trainInput; alternatives(selection, :)];
                obj.trainOutput = [obj.trainOutput; outputs(selection, :)];
                alternatives(selection, :) = [];
                outputs(selection, :) = [];
                stdDiff(nAlternatives, :) = [];
                nAlternatives = nAlternatives - 1;
            end
        end
        
        function activeLearnUpdate(obj, alternatives, outputs, nTest, nChanges)         
            testPoints = obj.inputSampler(nTest);         
            [~, prevStd] = obj.predict(testPoints);
            [nAlternatives, ~] = size(alternatives);
            [nTest, ~] = size(testPoints);
            stdDiff = zeros([nAlternatives, nTest]);
            tempPrediction = @(alt) obj.predictUpdated(testPoints, alt, 1);
            for iChange = 1:nChanges
                for iAlternative = 1:nAlternatives
                    [~, tempStd] = tempPrediction(alternatives(iAlternative, :));
                    stdDiff(iAlternative, :) = prevStd - tempStd;
                end
                stdDiffErr = std(stdDiff, 1, 2);
                stdDiff2 = mean(stdDiff, 2); 
                [~, selection] = max(stdDiff2);
                
                obj.trainInput = [obj.trainInput; alternatives(selection, :)];
                obj.trainOutput = [obj.trainOutput; outputs(selection, :)];
                alternatives(selection, :) = [];
                outputs(selection, :) = [];
                stdDiff(nAlternatives, :) = [];
                nAlternatives = nAlternatives - 1;
            end
        end   
        
        function obj = addNewPoints(obj, nTest, nChanges)
            alternatives = obj.libraryInputs;
            outputs = obj.libraryOutputs;
            for iChange = 1:nChanges
                [selections, ~, ~] = obj.activeLearn(alternatives, nTest);
                obj = obj.updateTrainingData(alternatives(selections, :), outputs(selections, :));
                alternatives(selections, :) = [];
                outputs(selections, :) = [];
            end
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
        function [selections, stdDiff, stdDiffErr] = activeUnlearn(obj, nTest, nChanges)
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
            [~, selections] = mink(stdDiff, nChanges);
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
       
       function newObj = copy(obj)
           % makes a shallow copy of the dynamic GP object.
           newObj = dynamicGP(obj.gpFuncHandle, obj.trainInput,...
            obj.trainOutput, obj.inputSampler, obj.libraryInputs,...
            obj.libraryOutputs);
       end
    end
end

