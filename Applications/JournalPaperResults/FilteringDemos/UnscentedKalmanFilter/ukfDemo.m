clear all
close all
clc
slope = 0.5;
additive = 5;
dataSlope = 3;

set(0, 'defaultFigureUnits', 'normalized')
set(0, 'defaultFigurePosition', [-0.4941    0.0312    0.4863    0.6111])

discrete = linspace(-35, 60, 1000);

initialState = 7; % 3 or 7
initialUncertainty = 5;

trueState = 7;
mappedTrueState = stateEvo(trueState, slope, additive);

figure(2)
subplot(2,2,4)
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
xlabel('state')
ylabel('pdf')
xlim([-35, 60])

subplot(2,2,2)
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
xlim([-35, 60])
ylim([-35,60])
xlabel('State')
ylabel('Evolved State')
hold on

disp('a')

stdDev = initialUncertainty;
samples = [initialState - 2*stdDev, initialState, initialState + 2*stdDev];
subplot(2,2,4)
scatter(samples, zeros(size(samples)), 120, 'black', 'o', 'filled')

disp('a')

for idx = 1:3
    %samples(idx) = normrnd(initialState, initialUncertainty);
    
    subplot(2,2,4)
    cla
    plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
    hold on
    scatter(samples, zeros(size(samples)), 120, 'black', 'o', 'filled')
    

    subplot(2,2,2)
    plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
    hold on
    scatter(samples(1:idx), -35*ones(size(samples(1:idx))), 'black', 'o')
    hold on
    plot([samples(idx), samples(idx)], [-35, stateEvo(samples(idx), slope, additive)], 'black', 'linewidth', 0.5)
    hold on
    scatter(-35*ones(size(samples(1:idx))), stateEvo(samples(1:idx), slope, additive), 120, 'black', 'x')
    hold on
    plot([-35, samples(idx)], [stateEvo(samples(idx), slope, additive), stateEvo(samples(idx), slope, additive)], 'black', 'linewidth', 0.5)
    
end
subplot(2,2,1)
scatter(stateEvo(samples, slope, additive), zeros(size(samples)), 120, 'black', 'o', 'filled');
hold on
xlim([-35,60])
ylim([0, 0.2])
legend("UKF Samples")
camroll(90)

disp('a')

subplot(2,2,1)
xlabel('state')
ylabel('pdf')
hold on
newSamples = normrnd(initialState, initialUncertainty^0.5, 1, 5000);
histogram(stateEvo(newSamples, slope, additive), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
hold on 
legend("Monte Carlo Approx", "UKF Samples")



function state = stateEvo(state, slope, add)
    state = slope*state + exp(0.2*state) + add - 2;
end

function data = dataGen(state)
    data = 3 * state;
end