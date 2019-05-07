function [samples, indexSet] = constructLibrary(largeTargetLogPdf, ...
    smallTargetLogPdf, ...
    sobolOptions,...
    bounds)

%% Input description
% largeTargetLogPdf                                                         Must take in nSamples x nDim matrix.  Return nSamples x 1 set of log density evaluations
% smallTargetLogPdf                                                         Same as above
% large target refers to the most uninformative distribution a filter can see
% small target refers to the most informative prior distribution a filter can see

%% sobolOptions object description
% contains parameters for determining how to construct offline library
% sobolOptions.terminationCriterion                                         Type of termination criterion, string, see minimumESS, or maxCardinality
% sobolOptions.minimumESS                                                   termitation criterion: minimum, scalar
% sobolOptions.maxCardinality                                               termination criterion: maximum number of samples in the library, scalar
% sobolOptions.thresholds                                                   set of density thresholds, vector of scalars
% sobolOptions.coolingStyle                                                 Metric to count for cooling times, string, 
% sobolOptions.coolingTimes                                                 count between threshold values, scalar

disp('starting bounds')
disp(bounds)
[~, numDim] = size(bounds);

[samples, indexSet] = validLibraryTargets(sobolOptions.thresholds,...
    numDim, ...
    sobolOptions.maxCardinality,...
    largeTargetLogPdf, bounds);


    function [library, indices] = validLibraryTargets(threshold, nDim,...
            maxCardinality, largeTargetLogPdf, bounds)
        ss = sobolset(nDim);
        iter = 1;
        nSamples = 0;
        while nSamples <= maxCardinality
            fprintf('Num samples %d \n', nSamples)
            if nSamples == 0
                propSamples = ss(1:maxCardinality, :);
                propSamples = convertSobolToBounds(propSamples, bounds);
                logP = largeTargetLogPdf(propSamples);
                indices = find(logP >= threshold);
                nSamples = length(indices);
            else
                propSamples = ss(iter*maxCardinality:(iter+1)*maxCardinality, :);
                propSamples = convertSobolToBounds(propSamples, bounds);
                logP = largeTargetLogPdf(propSamples);
                tempIndices = find(logP >= threshold);
                tempIndices = tempIndices + maxCardinality*ones(size(tempIndices));
                indices = [indices; tempIndices];
                nSamples = length(indices);
            end
        end
        library = convertSobolToBounds(ss(indices,:), bounds);
    end

    function samples = convertSobolToBounds(samples, bounds)
        % bounds =  2 x nDim matrix, lower bound is first row, upper = 2nd
        [nSamp, nDim] = size(samples);
        blankOnes = ones([nSamp, 1]);
        for iDim = 1:nDim
            samples(:,iDim) = samples(:,iDim)*(...
                bounds(2,iDim) - bounds(1,iDim)...
                ) + ...
                bounds(1,iDim)*blankOnes;
        end
    end
end