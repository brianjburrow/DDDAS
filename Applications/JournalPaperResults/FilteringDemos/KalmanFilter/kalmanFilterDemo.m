slope = 0.5;
additive = 2;
dataSlope = 3;

discrete = linspace(-35, 40, 1000);

initialState = 5;
initialUncertainty = 10;

trueState = 10;
mappedTrueState = stateEvo(trueState, slope, additive);
figure(1)
subplot(2,2,4)
plot(discrete, normpdf(discrete, initialState, initialUncertainty), 'black', 'linewidth', 2.0)
hold on
scatter(trueState, 0, 'red', 'x')
hold on
scatter(initialState, 0, 'black', 'o')
xlabel('state')
ylabel('pdf')
xlim([-35, 40])

subplot(2,2,2)
plot(discrete, stateEvo(discrete, slope, additive), 'blue', 'linewidth', 2.0, 'linestyle', '-.')
xlim([-35, 40])
xlabel('State')
ylabel('Evolved State')

subplot(2,2,1)
plot(discrete, normpdf(discrete, initialState*slope + additive, initialUncertainty*(slope.^2)), 'black', 'linewidth', 2.0)
xlabel('state')
ylabel('pdf')
hold on
scatter(mappedTrueState, 0, 'red', 'x')
hold on
scatter(stateEvo(initialState, slope, additive), 0, 'black', 'o')
xlim([-35,40])
camroll(90)


function state = stateEvo(state, slope, add)
    state = slope * state + add;
end

function data = dataGen(state)
    data = 3 * state;
end