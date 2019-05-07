%% EKF_demo
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
ylabel('Probability Density')

disp('a')

plot(xxx, ...
    normpdf(xxx, mean(initial_prior_samples), std(initial_prior_samples)),...
    'black', 'LineWidth', 2.0);

disp('a')

subplot(3,2,4)
poin = mean(initial_prior_samples);
y = cappe_ukfDataGen(poin);
slop = forward_slope(poin);
line = @(xx) slop*(xx - poin) + y;
plot(xxx, line(xxx), 'blue', 'LineWidth', 2.0)

disp('a')

subplot(3,2,2)

mu_x = mean(initial_prior_samples);
var_x = var(initial_prior_samples);
mu_y = cappe_ukfDataGen(mu_x);
var_y = var_x*slop^2;


% EKF Update
ekf_meas = line(initial_prior_samples);
d_til = mu_x;
e_til = mu_y;

data = [initial_prior_samples, prior_meas];
cross_cov = [var_x, 0; 0, var_y];

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

scatter(initial_prior_samples,...
        cappe_ukfDataGen_noisy(initial_prior_samples, 1), 30,...
        'filled')

    
disp('a')
    

scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')

subplot(3,2,1)

contour(X, Y, pdfVals_update)

disp('a')

scatter(mcmc_samples, mcmc_meas, 30,...
        'filled')

    
disp('a')
    
scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')


subplot(3,2,3)


plot(xxx, normpdf(xxx, d_t, P_t^0.5),...
    'black', 'LineWidth', 2.0)

disp('a')

histogram(mcmc_samples, 'Normalization', 'pdf')


disp('a')

subplot(3,2,5)
poin = mean(mcmc_samples);
y = cappe_ukfStateFunc(poin, 1);
slop = evolve_slope(poin);
line = @(xx) slop*(xx - poin) + y;
plot(xxx, line(xxx), 'blue', 'LineWidth', 2.0)

disp('a')

subplot(3,2,6)
cla
mu_x = mean(mcmc_samples);
var_x = var(mcmc_samples);
mu_y = cappe_ukfStateFunc(mu_x, 1);
var_y = var_x*slop^2;

plot(xxx, normpdf(xxx, mu_y, var_y^0.5), 'red', 'LineWidth', 2.0)

histogram(cappe_ukfStateFunc_noisy(mcmc_samples,1), 'Normalization', 'pdf')


function derivative = forward_slope(x)
    derivative = 2*0.05.*x;
end

function der = evolve_slope(x)
    top = x.^4 - 45.*x.^2 + 51;
    bottom = 2*(1 + x.^2).^2;
    der = top/bottom;
end