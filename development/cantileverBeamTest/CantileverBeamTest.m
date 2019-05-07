function [strains, M] = CantileverBeamTest(inputs)
    % this will be used to test how noise in the loading conditions will
    % affect the distribution of strain measurements that we get (i.e.,
    % how should we model the noise).
    
    % Inputs are [nSamples x 1001] matrix
    % columns 1:1000 are applied loads, the last column is the location
    % Set static load locations
    loadLocations = linspace(0, 1, 1000);
    L = 5000;
    E = 200000;
    I = 78125000;
    
    [nSamples, ~] = size(inputs);
    M = zeros([nSamples, 1]);
    for iSample = 1:nSamples
        loc = inputs(iSample, end);
        if (loc > L) || (loc < 0)
            error("location must be between 0 and 5000")
        end 
        indx = find(loadLocations > loc, 1);
        tempLoadLocs = loadLocations(indx:end);
        tempDist = tempLoadLocs - loc * ones(size(tempLoadLocs));
        temploads = inputs(iSample, indx:1000);
        M(iSample) = tempDist * temploads';
    end
    strains = (10^6)*M/(E*I);
end
        