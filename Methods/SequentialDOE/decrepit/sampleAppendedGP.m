function [alternativeGP, nAlternatives] = sampleAppendedGP(...
    GP, testSet, currentTrainingIn, currentTrainingOut,...
    newPoints, nSamples)
    %% sampleAppendedGP.mat
    % Author: Brian Burrows
    % The goal of this function is to support sequential experimental design
    % by takeing a function handle that makes predictions using a GPML
    % object.  We append the new datapoint, and sample the resulting
    % GP function nSamples times.
    
    % INPUT:
    % GP                                                                    GPML function handle, takes in three arguments: trainingSetIn, trainingSetOut, testSet
    % testSet                                                               nSamples x nDimensions array containing the total set of alternatives (or where you want to evaluate the GP)
    % currentTrainingSet                                                    Contains the array that we trained the GP on previously
    % newPoints                                                             nAlternatives x nDimensions array containing the alternative that you are considering adding to the GP
    % nSamples                                                              number of Monte Carlo samples desired to compute average entropy over the testSet
    

    %% Use Old GP to Determine Marginal distributions at each alternative 
    % used for sampling the next GP 
    [oldMu, oldStd] = GP(currentTrainingIn, currentTrainingOut, newPoints); % Make predictions using previous GP, used for entropy calculation
    
    [nAlternatives, ~] = size(newPoints);
    newMu = zeros(nAlternatives, nSamples);
    newStd = zeros(nAlternatives, nSamples);
    [nTestPoints, ~] = size(testSet);
    alternativeGP(nAlternatives).stdSet = zeros([nTestPoints, nSamples]);
    alternativeGP(nAlternatives).muSet = zeros([nTestPoints, nSamples]);
    for iAlternative = 1:nAlternatives
        if mod(iAlternative, 100) == 0
            fprintf("Analysing Alternative %d - %d \n", iAlternative, (iAlternative + 100))
        elseif iAlternative == 1
            fprintf("Analysing Alternative %d - %d \n", iAlternative, (iAlternative + 100))
        end
        possibleOutputs = normrnd(... 
            oldMu(iAlternative),...
            oldStd(iAlternative),...
            [nSamples,1]);
        
        newGP = @(zz) GP(...
                [currentTrainingIn; newPoints(iAlternative,:)],...
                [currentTrainingOut; zz],...
                testSet);   
            
        alternativeGP(iAlternative).muSet  = zeros([nTestPoints, nSamples]);
        alternativeGP(iAlternative).stdSet = zeros([nTestPoints, nSamples]);
        
        for iSample = 1:nSamples
            [alternativeGP(iAlternative).muSet(:, iSample), ...
                alternativeGP(iAlternative).stdSet(:, iSample)] = newGP(...
                possibleOutputs(iSample)...
                );
            alternativeGP(iAlternative).outputSet = currentTrainingOut; 
        end
    end
end