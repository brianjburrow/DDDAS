clear all
clc
close all

a = load('sample_set_1.mat');
b = load('sample_set_2.mat');

dam = load('Stored_data/lhs_4_param_6000_filter_states.txt');
dam = dam(:,1:4);

mcmc_samples_1 = a.mcmc_samples;
mcmc_samples_2 = b.mcmc_samples;

mu = zeros(99, 4);
errs = zeros(99, 4);

CAPS = transpose(load("capabilityTestData.npy"));
load('Surrogates/capfunc.mat');

for idx = 1:99
    bb = mcmc_samples_2(idx).val(10000:50:end,:);
    mu(idx, 1:4) = mean(bb);
    covs(idx, 1:4) = std(bb);
    mcmc_samples(idx).val = bb;
    [caps(idx).val, std_caps(idx).val] = predict(capFunc, bb);
    mu_caps(idx) = mean(caps(idx).val);
end
    
c = @colors;
name = 'maroon (x11)';

for dmx = 1:4
    subplot(2, 2, dmx)
    title("Bayesian Solution")
    hold on
    plot(mu(:,dmx),'color',c(name), 'LineWidth', 2.0)
    hold on
    plot(dam(:,dmx), 'black', 'LineWidth', 2.0, 'LineStyle', '-.');
    hold on
    errorbar(mu(:,dmx), 3*covs(:,dmx),'color',c(name), 'LineWidth', 0.5)
    hold on
    ylim([0,1])
    hold on
    xlabel('Time')
    hold on
    ylab = sprintf("State %d", dmx);
    ylabel(ylab)
end

for idx = 1:length(mu_caps)
    sig_caps(idx) = sum(std_caps(idx).val.^2) + sum(caps(idx).val.^2);
    sig_caps(idx) = (sig_caps(idx)/length(caps(idx).val) -...
        mean(caps(idx).val)^2)^0.5;
end


figure(2)
title("Bayesian Solution")
hold on
plot(mu_caps, 'color',c(name), 'LineWidth', 2.0);
hold on
errorbar(mu_caps, 3*sig_caps,'color',c(name), 'LineWidth', 0.5)
hold on
scatter(1:length(CAPS(2:end)), CAPS(2:end), 'k', 'filled')
legend({'Bayesian Estimate', 'Bayesian Uncertainty', 'True Capability',}, 'FontSize', 16)
xlabel("Time")
ylabel("Maximum Safe Load Factor")
ylim([4.75, 5.1])


