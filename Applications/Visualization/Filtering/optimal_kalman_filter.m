%% optimal_kalman_filter
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
cla
scatter(initialState, 0, 150, 'red', 'filled')
scatter(mean(initial_prior_samples), 0, 35, 'black')
disp('a')

histogram(initial_prior_samples, 'Normalization', 'pdf', 'FaceColor', 'blue')
ylabel('Probability Density')

disp('a')

subplot(3,2,2)
scatter(initial_prior_samples,...
        noisy_meas, 30,...
        'filled', 'blue')
   
disp('a')

subplot(3,2,2)
contour(X, Y, pdfVals)
disp('a')
    

scatter(initialState, measurement, 150, 'red', 'filled')

    
disp('a')

subplot(3,2,1)

contour(X, Y, pdfVals_update)
    
disp('a')
    
scatter(initialState, measurement, 150, 'red', 'filled')

disp('a')


subplot(3,2,3)
histogram(mcmc_samples, 'Normalization', 'pdf')
plot(xxx, normpdf(xxx, d_t, P_t^0.5), 'black', 'LineWidth', 2.0)


plot(xxx, normpdf(xxx, mean(mcmc_samples), std(mcmc_samples)), 'blue')

disp('a')

post_samples = normrnd(d_t, P_t^0.5, [1, 1000]);
histogram(post_samples, 'Normalization', 'pdf');

disp('a')

mapped_samples = cappe_ukfStateFunc_noisy(post_samples, 1);
subplot(3,2,6)
cla
histogram(cappe_ukfStateFunc_noisy(mcmc_samples,1), 'Normalization', 'pdf')
histogram(mapped_samples, 'Normalization', 'pdf')


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
