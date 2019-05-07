function testSampleAppendedGP(alternativeGPSet, testSet, nAlternatives, nSamples)
    [~, nDimensions] = size(testSet);
    if nDimensions ~= 1
        error("This function only works on 1-dimensional input spaces")
    end
    nFig = get(gcf,'Number');
    figure(nFig + 1)
    for iAlternative = 1:10:nAlternatives
        fprintf("Plotting Animation for Alternative %d \n", iAlternative)
        for iSample = 1:nSamples
            mu = alternativeGPSet(iAlternative).muSet(:, iSample);
            std = alternativeGPSet(iAlternative).stdSet(:, iSample);
            niceErrorbar(testSet, mu, std)
            drawnow
            pause(0.05)
            cla
        end
    end
end