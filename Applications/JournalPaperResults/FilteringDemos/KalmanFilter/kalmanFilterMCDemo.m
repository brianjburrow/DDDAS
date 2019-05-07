clear all
close all
clc
slope = 0.5;
additive = 5;
dataSlope = 3;

discrete = linspace(-35, 40, 1000);

initialState = 5;
initialUncertainty = 5;

trueState = 6;

figure(2)
subplot(2,2,4)
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
xlabel('state')
ylabel('pdf')
xlim([-35, 40])

subplot(2,2,2)
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
xlim([-35, 40])
ylim([-35,40])
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
    xlabel('state')
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