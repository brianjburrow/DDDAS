%% L2O_filter_demo
% NOTE: Bayes_filtering_demo.m must be run first in order for this to
% run correctly.

% This code was used to generate a demonstration plot for a presentation.
% There are a series of these demos.  They are for demonstrating concepts,
% and I do not necessarily use the correct method of computing some of
% these quantities.  In particular, the update stage is typically incorrect
% in each of the Kalman filtering demos (the cross-correlation dissappears
% when it should not), and the reason is that the Kalman filtering
% equations do not give you the joint posterior representation, it is
% merely useful for showing conceptually how the forward model factors in.


%% Prepare Figure

subplot(3,3,1)
title("Posterior Joint")
hold on
xlim([-50, 50])
ylim([-40, 100])
xlabel('State')
ylabel('Measurement')
subplot(3,3,2)
title("Prior Joint")
hold on
% Draw the error ellipse
xlabel('State')
ylabel('Measurement')

xlim([-50, 50])
ylim([-40, 100])
subplot(3,3,4)
title("State Estimate (Marginal)")
hold on
xlim([-50, 50])
xlabel('State')
ylabel('Probability Density')
subplot(3,3,5)
title("Forward Model")
hold on
xlim([-50, 50])
ylim([-40, 100])
xlabel('State')
ylabel('Measurement')
hold on
plot(xxx, yyy, 'black', 'LineWidth', 2.0)
subplot(3,3,7)
title("State Evolution Model")
hold on
plot(xxx, cappe_ukfStateFunc(xxx, 1), 'black', 'LineWidth', 2.0);
xlabel('State')
ylabel('State')
xlim([-50, 50])
subplot(3,3,8)
title("Predicted Prior")
xlabel('State')
hold on
xlim([-50, 50])
ylim([0, inf])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.98, 0.74])

subplot(3,3,3)
title("Weighted Joint")
xlim([-50, 50])
ylim([-40, 100])
hold on
subplot(3,3,6)
title("Offline Forward Evals")
hold on
olib_x = linspace(-50, 50, 50);
olib_y = cappe_ukfDataGen(olib_x);
scatter(olib_x, olib_y, 10, 'black', 'filled')
xlim([-50, 50])
ylim([-40, 100])
subplot(3,3,9)
title("Proposal CDF")
hold on
ecdf(olib_x)
xlim([-50, 50])
ylim([0, 1])


%% Update Figure

subplot(3,3,8)
cla
scatter(initialState, 0, 150, 'red', 'filled')
scatter(mean(initial_prior_samples), 0, 35, 'black')
histogram(initial_prior_samples, 'Normalization', 'pdf', 'FaceColor', 'blue')
ylabel('Probability Density')

disp('a')

cla
ecdf(initial_prior_samples);

disp('a')
    
subplot(3,3,9)
weights = plot_l2o_update(olib_x, initial_prior_samples');
weights = 1250*weights;
disp('a')

subplot(3,3,3)
scatter(olib_x, olib_y, weights, 'blue', 'filled')


subplot(3,3,2)
scatter(initial_prior_samples,...
        noisy_meas, 10,...
        'filled', 'blue')
   
disp('a')

subplot(3,3,3)
contour(X, Y, pdfVals)
    disp('a')
    

scatter(initialState, measurement, 150, 'red', 'filled')

    
disp('a')

subplot(3,3,1)

contour(X, Y, pdfVals_update)
    
disp('a')
    
scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')


subplot(3,3,4)
histogram(mcmc_samples, 'Normalization', 'pdf')
plot(xxx, normpdf(xxx, d_t, P_t^0.5), 'black', 'LineWidth', 2.0)

disp('a')

subplot(3,3,8)
cla
newPriorSamples = cappe_ukfStateFunc_noisy(mcmc_samples',1);
histogram(newPriorSamples, 'Normalization', 'pdf')
disp('a')
cla
ecdf(newPriorSamples')
disp('a')
subplot(3,3,9)
cla
ecdf(olib_x)
weights = plot_l2o_update(olib_x, newPriorSamples);
weights = 250*weights;
disp('a')

subplot(3,3,2)
cla
scatter(newPriorSamples, cappe_ukfDataGen_noisy(newPriorSamples), 3, 'filled')

subplot(3,3,3)
cla
scatter(olib_x, olib_y, weights, 'blue', 'filled')
