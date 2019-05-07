clear all
close all
clc

set(0, 'defaultFigureUnits', 'normalized')
set(0, 'defaultFigurePosition', [-0.4941    0.0312    0.4863    0.6111])
slope = 0.5;
additive = 5;
dataSlope = 3;
plotUB = 100;
measurementNoise = 10;

discrete = linspace(-35, plotUB, 1000);

initialState = -3; % 3 or 7
initialUncertainty = 10;

trueState = 14;
mappedTrueState = stateEvo(trueState, slope, additive);

figure(1)
hold on
subplot(2,2,4)
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
ylabel('pdf')
xlim([-35, plotUB])

disp('a')

subplot(2,2,2)
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
hold on
xlim([-35, plotUB])
ylim([-35,plotUB])
scatter(trueState, -35, 'blue', 'x')
hold on
scatter(initialState, -35, 'black', 'o')
ylabel('Evolved State')
hold on
drawnow

disp('a')

plot(discrete, linearState(discrete, initialState, slope, additive), 'blue', 'linewidth', 2.0)
hold on
drawnow

disp('a')
    

for idx = 1:100
    samples(idx) = normrnd(initialState, initialUncertainty);
    if idx == 100
        samples(idx) = 21;
    end
    subplot(2,2,4)
    cla
    plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
    hold on
    scatter(samples(idx), 0, 'black', 'o')
    drawnow

    subplot(2,2,2)
    cla()
    plot(discrete, linearState(discrete, initialState, slope, additive), 'blue', 'linewidth', 2.0)
    hold on
    plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
    hold on
    scatter(samples(1:idx), -35*ones(size(samples(1:idx))), 'black', 'o')
    hold on
    plot([samples(idx), samples(idx)], [-35, stateEvo(samples(idx), slope, additive)], 'black', 'linewidth', 0.5)
    hold on
    scatter(-35*ones(size(samples(1:idx))), stateEvo(samples(1:idx), slope, additive), 'green', 'x')
    hold on
    plot([-35, samples(idx)], [stateEvo(samples(idx), slope, additive), stateEvo(samples(idx), slope, additive)], 'green', 'linewidth', 0.5)
    hold on
    scatter(-35*ones(size(samples(1:idx))), linearState(samples(1:idx), initialState, slope, additive), 'black', 'o')
    hold on
    plot([-35, samples(idx)], ...
        [linearState(samples(idx), initialState, slope, additive),...
        linearState(samples(idx), initialState, slope, additive)],...
        'black', 'linewidth', 0.5, 'linestyle', '--')
 

    subplot(2,2,1)
    cla()
    ylabel('pdf')
    hold on
    histogram(stateEvo(samples(1:idx), slope, additive), 'Normalization', 'pdf', 'FaceColor', 'green', 'FaceAlpha', 0.2)
    hold on
    histogram(linearState(samples(1:idx), initialState, slope, additive), 'Normalization', 'pdf', 'FaceColor', 'black', 'FaceAlpha', 0.2)
  
    if idx == 1
        xlim([-35,plotUB])
        camroll(90)
    end
end

disp('a')

subplot(2,2,1)
plot(discrete, ...
    normpdf(discrete, mean(linearState(samples,  initialState, slope, additive)), std(linearState(samples, initialState,slope,additive))),...
    'black', 'linewidth', 2.0)
mtit('Predict Step')

newMean = mean(linearState(samples, initialState, slope, additive));
newStd = std(linearState(samples,initialState, slope,additive));

disp('a')

figure(2)
hold on
subplot(2,2,1)
plot(discrete, ...
    normpdf(discrete, mean(linearState(samples, initialState, slope, additive)), std(linearState(samples,initialState,slope,additive))),...
    'black', 'linewidth', 2.0)
hold on 
scatter(min(dataGen(discrete)), stateEvo(trueState, slope, additive), 'blue', 'x')
hold on
xlim([-35, plotUB])
camroll(90)
ylabel('PDF')
title('Prior distribution')

subplot(2,2,2)
xlabel("Measurements")
ylabel("State")
hold on
ylim([-35, plotUB])
hold on
xlim([dataGen(min(discrete)), dataGen(plotUB)])
hold on
plot(dataGen(discrete), discrete, 'blue', 'linewidth', 2.0,'linestyle','-.')
hold on
plot([dataGen(mappedTrueState), dataGen(mappedTrueState)], [mappedTrueState, min(dataGen(discrete))], 'blue', 'linewidth', 1.0,'linestyle','-')
hold on
plot([min(dataGen(discrete)), dataGen(mappedTrueState)], [mappedTrueState, mappedTrueState], 'blue', 'linewidth', 1.0,'linestyle','-')
hold on
scatter([min(dataGen(discrete)), min(dataGen(discrete))], [newMean, newMean], 'black', 'o')
hold on 
scatter(min(dataGen(discrete)), stateEvo(trueState, slope, additive), 'blue', 'x')
hold on
scatter(dataGen(mappedTrueState), min(dataGen(discrete)), 'blue', 'x')
hold on
plot([min(dataGen(discrete)), dataGen(newMean)], [newMean, newMean], 'black')
hold on
plot([dataGen(newMean), dataGen(newMean)],[newMean, -35], 'black');
hold on
scatter([dataGen(newMean), dataGen(newMean)],[-35, -35], 'black', 'x');
title('Likelihood')

subplot(2,2,4)
plot(dataGen(discrete),normpdf(dataGen(discrete), dataGen(newMean), newStd*3), 'black', 'linewidth', 2.0)
hold on
scatter(dataGen(mappedTrueState), 0, 'blue', 'x')
xlabel('Measurements')
xlim([-100, 300])

newMean = newMean;

mu = (1 / ((1/newStd^2) + (1/measurementNoise))) * (dataGen(newMean)/newStd^2 + dataGen(mappedTrueState)/measurementNoise);
mu = mu/3;
sigma = ((1/newStd^2) + (1/measurementNoise))^-1;

subplot(2,2,1)
plot(discrete, ...
    normpdf(discrete, mu, sigma^0.5),...
    'black', 'linewidth', 2.0, 'linestyle', '-.')


function state = linearState(state, setState, slope, add)
    y1 = stateEvo(setState + 0.000001, slope, add);
    y2 = stateEvo(setState - 0.000001, slope, add);
    
    slope = (y2 - y1)/(-2*0.000001);  
    
    intercept = y1 - slope*(setState + 0.000001);
   
    state = state * slope + intercept;
end
function state = stateEvo(state, slope, add)
    state = slope*state + exp(0.2*state) + add - 2;
end

function data = dataGen(state)
    data = 3 * state;
end

    