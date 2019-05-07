clear all
close all
clc

%% propagationDemo1
% This function is designed to show how Monte Carlo simulation can be used
% to propagate a distribution through a function.  We generate samples from
% the initial distribution, evaluate each sample through the function of
% interest.  Think y_i = f(x_i), where x_i are samples from our initial
% distribution p(x).  Then the set of samples {y_i for i in 1,...,k} forms
% a k sample approximation of our output distribution of interest, p(y).
% Also, Gaussian distributions remain Gaussian after a linear
% transformation (See propagationDemo2 for more info).
set(0, 'defaultFigureUnits', 'normalized')
set(0, 'defaultFigurePosition', [-0.4941    0.0312    0.4863    0.6111])

slope = 0.5;
additive = 5;
dataSlope = 3;

discrete = linspace(-35, 40, 1000);

initialState = 5;
initialUncertainty = 5;

trueState = 10;
mappedTrueState = stateEvo(trueState, slope, additive);

figure(2)
subplot(2,2,4)
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
xlabel('state')
ylabel('pdf')
xlim([-35, 40])

disp('a')

subplot(2,2,2)
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
xlim([-35, 40])
ylim([-35,40])
xlabel('State')
ylabel('Evolved State')
hold on
drawnow
idx = 1;
samples(idx) = normrnd(initialState, initialUncertainty);

disp('a')

subplot(2,2,4)
cla
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
scatter(samples(idx), 0, 'black', 'o')

disp('a')

subplot(2,2,2)
cla()
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
hold on
scatter(samples(1:idx), -35*ones(size(samples(1:idx))), 'black', 'o')
hold on
plot([samples(idx), samples(idx)], [-35, stateEvo(samples(idx), slope, additive)], 'black', 'linewidth', 0.5)
hold on
scatter(-35*ones(size(samples(1:idx))), stateEvo(samples(1:idx), slope, additive), 'red', 'x')
hold on
plot([-35, samples(idx)], [stateEvo(samples(idx), slope, additive), stateEvo(samples(idx), slope, additive)], 'red', 'linewidth', 0.5)

disp('a')


subplot(2,2,1)
cla()
ylabel('pdf')
hold on
histogram(stateEvo(samples(1:idx), slope, additive), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
    
for idx = 1:100
    samples(idx) = normrnd(initialState, initialUncertainty);
    
    subplot(2,2,4)
    cla
    plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
    hold on
    scatter(samples(idx), 0, 'black', 'o')
    
    drawnow

    subplot(2,2,2)
    cla()
    plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
    hold on
    scatter(samples(1:idx), -35*ones(size(samples(1:idx))), 'black', 'o')
    hold on
    plot([samples(idx), samples(idx)], [-35, stateEvo(samples(idx), slope, additive)], 'black', 'linewidth', 0.5)
    hold on
    scatter(-35*ones(size(samples(1:idx))), stateEvo(samples(1:idx), slope, additive), 'red', 'x')
    hold on
    plot([-35, samples(idx)], [stateEvo(samples(idx), slope, additive), stateEvo(samples(idx), slope, additive)], 'red', 'linewidth', 0.5)
 

    subplot(2,2,1)
    cla()
    ylabel('pdf')
    hold on
    histogram(stateEvo(samples(1:idx), slope, additive), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
  
    if idx == 1
        xlim([-35,40])
        camroll(90)
    end
end

function state = stateEvo(state, slope, add)
    state = slope * state + add;
end

function data = dataGen(state)
    data = 3 * state;
end