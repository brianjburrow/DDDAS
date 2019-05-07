clear all
clc
close all

polyOrders = 1;
dim = 1;  % Dimensionality of parameter space
lb = 1;
ub = 10;
numSlice = 50;
numSamples = 200;
%samps = genPosteriorSamples(4);


%% Run Script
XXX = linspace(lb, ub, numSlice*2);
xxx = linspace(lb, ub, numSlice);
samples = zeros(numSamples, numSlice);
figure(1)
subplot(4,1,1)
for idx = 1:length(xxx)
    mu(idx) = 10*xxx(idx);
    Sigma(idx) = 10*abs(sin(xxx(idx)));
end

errorbar(xxx, mu, 3*Sigma.^0.5)

x1r = fitrsvm(xxx', mu,'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true, 'OptimizeHyperparameters','auto');

x2r = fitrsvm(xxx', Sigma,'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true, 'OptimizeHyperparameters','auto');


figure(1)
subplot(4,1,2)
plot(XXX, predict(x1r, XXX'), 'red')
hold on
scatter(XXX, 10*XXX);
drawnow 
xlim([1, 10])
ylim([0, 100])
hold on
xlabel("Input")
ylabel("P_{tar}")

figure(1)
subplot(4,1,3)
plot(xxx, predict(x2r, xxx'), 'red')
hold on
scatter(xxx, 10*abs(sin(xxx)));

figure(1)
subplot(4,1,4)
errorbar(xxx, zeros(size(xxx)), 3*ones(size(xxx)), 'red')
hold on
