clear all
close all
clc
slope = 0.5;
additive = 5;
dataSlope = 3;

discrete = linspace(-35, 60, 1000);

initialState = 3; % 3 or 7
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
drawnow
disp('a')
for idx = 1:1000
    samples(idx) = normrnd(initialState, initialUncertainty);
    
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
    scatter(-35*ones(size(samples(1:idx))), stateEvo(samples(1:idx), slope, additive), 'red', 'x')
    hold on
    plot([-35, samples(idx)], [stateEvo(samples(idx), slope, additive), stateEvo(samples(idx), slope, additive)], 'red', 'linewidth', 0.5)
    hold on
    scatter(-35*ones(size(samples(1:idx))), linearState(samples(1:idx), initialState, slope, additive), 'green', 'o')
    hold on
    plot([-35, samples(idx)], [linearState(samples(idx), initialState, slope, additive), linearState(samples(idx), initialState, slope, additive)], 'green', 'linewidth', 0.5)
 

    subplot(2,2,1)
    cla()
    xlabel('state')
    ylabel('pdf')
    hold on
    histogram(stateEvo(samples(1:idx), slope, additive), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.2)
    hold on
    histogram(linearState(samples(1:idx), initialState, slope, additive), 'Normalization', 'pdf', 'FaceColor', 'green', 'FaceAlpha', 0.2)
  
    if idx == 1
        xlim([-35,60])
        camroll(90)
    end
end

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