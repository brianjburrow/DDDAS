close all
clear all
clc
%% Define variables
%% bayes_filtering_demo
% This code was used to generate a demonstration plot for a presentation.
% There are a series of these demos.  They are for demonstrating concepts,
% and I do not necessarily use the correct method of computing some of
% these quantities.  In particular, the update stage is typically incorrect
% in each of the Kalman filtering demos (the cross-correlation dissappears
% when it should not), and the reason is that the Kalman filtering
% equations do not give you the joint posterior representation, it is
% merely useful for showing conceptually how the forward model factors in.

total_time = 20;   % Total number of time steps
chainLength = 20000;
%% Define Filters
% Initial Conditions
initialState = 0;
processNoise = 10;
measurementNoise = 1;
measurement = cappe_ukfDataGen_noisy(initialState, 1);
measurement = cappe_ukfDataGen_noisy(initialState, 1) + 10;
%measurement = cappe_ukfDataGen_noisy(initialState, 1) + 50;

% define unscented kalman filter
ukf = unscentedKalmanFilter( ...
@cappe_ukfStateFunc,...
@cappe_ukfDataGen, ...
initialState, ...
'HasAdditiveMeasurementNoise',1);
ukf.ProcessNoise = processNoise;
ukf.MeasurementNoise = measurementNoise;
ukf.Alpha = 1;
ukf.Kappa = 3;
ukf.Beta = 8;

% define extended kalman filter
%% Initialize EKF
ekf = extendedKalmanFilter( ...
@cappe_ukfStateFunc,...
@cappe_ukfDataGen, ...
initialState, ...
'HasAdditiveMeasurementNoise',1);
ekf.ProcessNoise = processNoise;
ekf.MeasurementNoise = measurementNoise;

initial_prior_samples = [normrnd(initialState+5, processNoise, [1000, 1])];
xxx = linspace(-60, 60, 1000);
yyy = cappe_ukfDataGen(xxx, 1);

prior_meas = cappe_ukfDataGen(initial_prior_samples, 1);

mu_x = mean(initial_prior_samples);
mu_y = mean(prior_meas);


data = [initial_prior_samples, prior_meas];
cross_cov = cov(data);

[X, Y] = ndgrid(-60:60, -120:2:120);
pdfVals = zeros(size(X));
pdfFunc = @(x,y) mvnpdf([x,y],[mu_x, mu_y], cross_cov);
for i=1:size(X,1)
    for j=1:size(X,2)
        pdfVals(i,j) = pdfFunc(X(i,j),Y(i,j));
    end
end

% Monte Carlo Update
P_e = cov(prior_meas);
P_til = cov(initial_prior_samples);
d_til = mu_x;
e_til = mu_y;
P_de = 0;
for idx = 1:length(initial_prior_samples)
    P_de = P_de + (initial_prior_samples(idx) - d_til) *...
        (prior_meas(idx) - e_til);
end
P_de = P_de/(length(initial_prior_samples));
Kt = P_de/P_e;
P_t = P_til - Kt*P_e*Kt';
d_t = d_til + Kt*(measurement - e_til);

P_meas = (1/cross_cov(2,2) + 1/measurementNoise)^-1;
mu_meas = (1/cross_cov(2,2) + 1/measurementNoise)*(mu_y/cross_cov(2,2) + ...
    measurement/measurementNoise);



cross_cov_update = zeros([2,2]);
cross_cov_update(1,1) = P_t;
cross_cov_update(2,2) = P_meas;

[X, Y] = ndgrid(-60:60, -120:2:120);
pdfVals_update = zeros(size(X));
pdfFunc_update = @(x,y) mvnpdf([x,y],[d_t, mu_meas], cross_cov_update);
for i=1:size(X,1)
    for j=1:size(X,2)
        pdfVals_update(i,j) = pdfFunc_update(X(i,j),Y(i,j));
    end
end
%% MCMC Estimate

mcmc_samples = zeros([chainLength, 1]);
targetPdf = @(xx) logmvnpdf(cappe_ukfDataGen(xx), measurement, 1) + logmvnpdf(xx, mu_x, 10);

C = 10;
mcmc_samples(1) = initialState;
for idx = 2:chainLength
    proposalState = normrnd(mcmc_samples(idx-1, 1), C);
    propPdf = @(xx, yy, zz) logmvnpdf(xx, mcmc_samples(idx-1, 1), zz);
    [mcmc_samples(idx, 1), accept, C] = metropolisHastingsIter(...
        mcmc_samples(idx-1, 1), ...
        proposalState,...
        propPdf, ...
        targetPdf, ...
        idx, ...
        C, ...
        mcmc_samples(1:idx, 1));
end
mcmc_samples = mcmc_samples(2000:20:end,1);
mcmc_meas = cappe_ukfDataGen_noisy(mcmc_samples, 1);

disp('a')


%% Prepare Figure

    subplot(3,2,1)
    title("Posterior Joint")
    hold on
    xlim([-50, 50])
    ylim([-40, 100])
    xlabel('State')
    ylabel('Measurement')
    subplot(3,2,2)
    title("Prior Joint")
    hold on
    % Draw the error ellipse
    xlabel('State')
    ylabel('Measurement')
    
    xlim([-50, 50])
    ylim([-40, 100])
    subplot(3,2,3)
    title("State Estimate (Marginal)")
    hold on
    xlim([-50, 50])
    xlabel('State')
    ylabel('Probability Density')
    subplot(3,2,4)
    title("Forward Model")
    hold on
    xlim([-50, 50])
    ylim([-40, 100])
    xlabel('State')
    ylabel('Measurement')
    hold on
    plot(xxx, yyy, 'black', 'LineWidth', 2.0)
    subplot(3,2,5)
    title("State Evolution Model")
    hold on
    plot(xxx, cappe_ukfStateFunc(xxx, 1), 'black', 'LineWidth', 2.0);
    xlabel('State')
    ylabel('State')
    xlim([-50, 50])
    subplot(3,2,6)
    title("Predicted Prior")
    xlabel('State')
    hold on
    xlim([-50, 50])
    ylim([0, inf])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.98, 0.74])

%% Update Figure

subplot(3,2,6)
scatter(initialState, 0, 150, 'red', 'filled')

disp('a')

histogram(initial_prior_samples, 'Normalization', 'pdf', 'FaceColor', 'blue')
ylabel('Probability Density')
disp('a')

subplot(3,2,2)

scatter(initialState, measurement, 150, 'red', 'filled')
disp('a')

noisy_meas = cappe_ukfDataGen_noisy(initial_prior_samples, 1);
scatter(initial_prior_samples,...
        noisy_meas, 30,...
        'filled', 'blue')
scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')

subplot(3,2,1)

scatter(mcmc_samples, mcmc_meas, 30,...
        'filled', 'blue')
    
disp('a')
    
scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')


subplot(3,2,3)
histogram(mcmc_samples, 'Normalization', 'pdf', 'FaceColor', 'blue')

plot(xxx, normpdf(xxx, mean(mcmc_samples), std(mcmc_samples)))

disp('a')

subplot(3,2,6)
cla
newPriorSamples = cappe_ukfStateFunc_noisy(mcmc_samples',1);
histogram(newPriorSamples, 'Normalization', 'pdf', 'FaceColor', 'blue')

disp('a')

subplot(3,2,2)
cla
scatter(newPriorSamples, cappe_ukfDataGen_noisy(newPriorSamples), 3, 'blue','filled')
% figure(2)
% scatterhist(initial_prior_samples,...
%         cappe_ukfDataGen_noisy(initial_prior_samples, 1))
% hold on
% xlim([-100, 100])
% ylim([0, 120])
% mu_x = mean(initialState);
% mu_y = mean(prior_meas);
% sig_x = std(initial_prior_samples);
% 
% figure(3)
% plot(xxx, normpdf(xxx, mu_x, cross_cov(1,1)^0.5))
% hold on
% histogram(initial_prior_samples, 'Normalization', 'pdf')
% xlim([-100, 100])
% 
% figure(4)
% plot(linspace(0, 120, 1000),...
%     normpdf(linspace(0, 120, 1000), mu_y, cross_cov(2,2)^0.5))
% hold on
% histogram(prior_meas, 'Normalization', 'pdf')
% xlim([0, 120])
% set(gca,'YDir','90');
% 
