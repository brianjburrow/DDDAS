function plot_mcmc_sampleset_caps(samples, truths)

    numChains = length(samples(1,:));
    numSamples = length(samples(:,1));
    disp('where am i')
    for idx = 1:numChains
        disp(idx)
        figure(idx)
        subplot(2,1,1)
        plot(samples(:,idx))
        hold on
        plot([0, numSamples], [truths(idx, 1), truths(idx, 1)], 'black', 'linewidth', 2.0)
        ylim([4.25,5.1])
        xlim([0, numSamples])

        
        subplot(2,1,2)
        histogram(samples(:,idx), 'normalization', 'probability')
        hold on
        plot([truths(idx, 1), truths(idx, 1)], [0, 1.], 'red', 'linewidth', 2.0, 'linestyle', '-.')
        ylim([0,0.01])
        xlim([4.25,5.1])

    end
end