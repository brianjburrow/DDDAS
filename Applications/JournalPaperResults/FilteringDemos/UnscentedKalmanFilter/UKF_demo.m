

%% UKF_demo
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
subplot(3,2,1)
title("Posterior Joint")
hold on
xlim([-50, 50])
ylim([-40, 160])
xlabel('State')
ylabel('Measurement')
subplot(3,2,2)
title("Prior Joint")
hold on
% Draw the error ellipse
xlabel('State')
ylabel('Measurement')

xlim([-50, 50])
ylim([-40, 160])
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
ylim([-40, 160])
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
ylim([-50, 50])
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

histogram(initial_prior_samples, 'Normalization', 'pdf')
disp('a')

plot(xxx, ...
    normpdf(xxx, mean(initial_prior_samples), std(initial_prior_samples)),...
    'black', 'LineWidth', 2.0);
disp('a')
sig_points = [-5.5, 5, 15.5];
histogram(sig_points, 30, 'Normalization', 'pdf')
disp('a')

ylabel('Probability Density')

disp('a')

subplot(3,2,2)

mu_x = mean(initial_prior_samples);
y_samps = cappe_ukfDataGen_noisy(sig_points);
mu_y = mean(y_samps);

scatter(sig_points, y_samps, 'black', 'filled');

disp('a')

data = [sig_points; y_samps];

cross_cov = cov(data');

[X, Y] = ndgrid(-60:60, -120:2:120);
pdfVals = zeros(size(X));
pdfFunc = @(x,y) mvnpdf([x,y],[mu_x, mu_y], cross_cov);
for i=1:size(X,1)
    for j=1:size(X,2)
        pdfVals(i,j) = pdfFunc(X(i,j),Y(i,j));
    end
end

contour(X, Y, pdfVals)

disp('a')  
    

scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')

subplot(3,2,1)


    
contour(X, Y, pdfVals_update)
    

disp('a')


subplot(3,2,3)


plot(xxx, normpdf(xxx, d_t, P_t^0.5), 'black', 'LineWidth', 2.0);

disp('a')

sigPointUp = [d_t - 2*P_t^0.5, d_t, d_t + 2*P_t^0.5];
histogram(sigPointUp, 35, 'Normalization', 'pdf')

disp('a')

subplot(3,2,6)
cla
sigPointFin = cappe_ukfStateFunc_noisy(sigPointUp, 1);

histogram(sigPointFin, 35, 'Normalization', 'pdf')

disp('a')
plot(xxx, normpdf(xxx, mean(sigPointFin), std(sigPointFin)), 'black', 'LineWidth', 2.0)

disp('a')

histogram(cappe_ukfStateFunc_noisy(mcmc_samples,1), 'Normalization', 'pdf')


function derivative = forward_slope(x)
    derivative = 2*0.05.*x;
end

function der = evolve_slope(x)
    top = x.^4 - 45.*x.^2 + 51;
    bottom = 2*(1 + x.^2).^2;
    der = top/bottom;
end