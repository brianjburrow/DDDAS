function plot_mcmc_sampleset(samples, truths)

    numChains = length(samples(:,1,1));
    numSamples = length(samples(1,:,1));
    numDimensions = length(samples(1,1,:));
    
    for idx = 1:numChains
        figure(idx)
        subplot(2,4,1)
        plot(samples(idx,:,1))
        hold on
        plot([0, numSamples], [truths(idx, 1), truths(idx, 1)], 'black', 'linewidth', 2.0)
        ylim([0,1])
        xlim([0, numSamples])
        subplot(2,4,2)
        plot(samples(idx,:,2))
        hold on
        plot([0, numSamples], [truths(idx, 2), truths(idx, 2)], 'black', 'linewidth', 2.0)
        ylim([0,1])
        xlim([0, numSamples])
        subplot(2,4,3)
        plot(samples(idx,:,3))
        hold on
        plot([0, numSamples], [truths(idx, 3), truths(idx, 3)], 'black', 'linewidth', 2.0)
        ylim([0,1])
        xlim([0, numSamples])
        subplot(2,4,4)
        plot(samples(idx,:,4))
        hold on
        plot([0, numSamples], [truths(idx, 4), truths(idx, 4)], 'black', 'linewidth', 2.0)
        ylim([0,1])
        xlim([0, numSamples])
        
        subplot(2,4,5)
        histogram(samples(idx,:,1), 'normalization', 'probability')
        hold on
        plot([truths(idx, 1), truths(idx, 1)], [0, 1.], 'red', 'linewidth', 2.0, 'linestyle', '-.')
        ylim([0,0.2])
        xlim([0,1])
        subplot(2,4,6)
        histogram(samples(idx,:,2), 'normalization', 'probability')
        hold on
        plot([truths(idx, 2), truths(idx, 2)], [0, 1.], 'red', 'linewidth', 2.0, 'linestyle', '-.')
        ylim([0,0.2])
        xlim([0,1])
        subplot(2,4,7)
        histogram(samples(idx,:,3), 'normalization', 'probability')
        hold on
        plot([truths(idx, 3), truths(idx, 3)], [0, 1.], 'red', 'linewidth', 2.0, 'linestyle', '-.')
        ylim([0,0.2])
        xlim([0,1])
        subplot(2,4,8)
        histogram(samples(idx,:,4), 'normalization', 'probability')
        hold on
        plot([truths(idx, 4), truths(idx, 4)], [0, 1.], 'red', 'linewidth', 2.0, 'linestyle', '-.')
        ylim([0,0.2])
        xlim([0,1])
    end
end