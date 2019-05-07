clear all
close all
clc
xx = linspace(-20, 20, 1000);
samples = normrnd(5,5,10000, 1);
bbSamples = nonlinearFunc(samples, 1);
bbSamples = cappe_ukfStateFunc(samples,1);
d = -inf;
e = inf;
n = 1;
lambda = 2;
sigP = SPICUKF_sigma_ICUT(5, 5^2, d, e, n, lambda);
ccSamples = nonlinearFunc(sigP, 1);
ccSamples = cappe_ukfStateFunc(sigP, 1);
ccSamples2 = cappe_ukfStateFunc(sigP, 1);

%sigP2 = SPICUKF_sigma_ICUT(5, 5^2, d, e, n, 0);
%ccSamples2 = cappe_ukfStateFunc(sigP2, 1);


subplot(3, 2, 1);
histogram(sigP, 35, 'FaceColor', 'black')
hold on
%histogram(sigP2, 35, 'FaceColor', 'blue')
xlabel('State', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)
ylim([0,3])
xlim([-10, 30])
title('Initial Prior', 'FontSize', 16)

subplot(3, 2, 2);
histogram(samples, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
hold on
xlabel('State', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16)
title('Initial Prior', 'FontSize', 16)
xlim([-10, 30])

disp('a')


subplot(3, 2, 3);
histogram(ccSamples, 35, 'FaceColor', 'black')
hold on
%histogram(ccSamples2, 35, 'FaceColor', 'blue')
hold on
ylim([0,3])
xlabel('State', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)
title('Predicted Prior', 'FontSize', 16)
xlim([-20, 20])

subplot(3, 2, 4);
histogram(bbSamples, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
hold on
hold on
xlabel('State', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16)
xlim([-20, 20])
title('Predicted Prior', 'FontSize', 16)

disp('a')


subplot(3, 2, 4);
plot(xx, normpdf(xx, mean(ccSamples), var(ccSamples)^0.5), 'black', 'LineWidth', 2.0)
hold on
xlabel('State', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16)
title('Initial Prior', 'FontSize', 16)

disp('a')

subplot(3, 2, 5);
histogram(ccSamples.^2, 35, 'FaceColor', 'black')
hold on
hold on
xlabel('State', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16)
xlim([0, 300])
ylim([0, 3])
title('Predicted Measurement', 'FontSize', 16)


subplot(3,2,6)
hh = histogram(bbSamples.^2, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2);
hold on
mu = mean(bbSamples);
va = var(bbSamples);
yy = normpdf(xx, mu, va^0.5);
xlim([0,300])
title('Predicted Measurements', 'FontSize', 16)


xlabel('State', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16)
%legend('UKF Histogram: \lambda_1', 'UKF Histogram: \lambda_2', 'MC Histogram', 'MC Mean', 'UKF Mean: \lambda_1', 'UKF Mean: \lambda_2')

disp('a')

subplot(3,2,6)
plot(xx.^2, normpdf(xx.^2, mean(ccSamples.^2), var(ccSamples.^2)^0.5), 'black', 'LineWidth', 2.0);
xlim([0,300])
title('Predicted Measurements', 'FontSize', 16)


function yy = nonlinearFunc(xx, dummy)
    yy = sin(xx) + 0.24*xx.^2;
end
