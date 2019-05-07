function sortIndices = determine_map_ordering(samples)
    [numSamples, numDim] = size(samples);
    Gamma = zeros([numDim, numDim]);
    theta_bar = mean(samples);
    for idx = 1:numDim
        for dmx = 1:numDim
            for tmz = 1:numSamples
                Gamma(idx, dmx) = Gamma(idx, dmx) + ...
                    (samples(tmz,idx) - theta_bar(idx)) *...
                    (samples(tmz, dmx) - theta_bar(dmx))^2;
            end
        end
    end
    
    Gamma = Gamma./numSamples;
    
    sortIndices = 1:numDim;
    
    for idx = 2:numDim
        j = idx;
        while (j > 1) && (Gamma(j, j-1) < Gamma(j-1, j))
            s1 = sortIndices(j-1);
            s2 = sortIndices(j);
            sortIndices(j-1) = s2;
            sortIndices(j) = s1;
        end
    end
end