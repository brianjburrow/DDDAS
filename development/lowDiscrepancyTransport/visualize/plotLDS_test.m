clear all
clc
close all

set(0, 'defaultFigureUnits', 'normalized')
set(0, 'defaultFigurePosition', [-0.4941    0.0312    0.8863    0.8111])


ss = sobolset(1);
numMax = 2000;
yy = funct(linspace(0,1,100000));

truth = sum(yy)/length(yy);

figure(1)
subplot(2,3,1)
xlim([0,1])
hold on
ylim([0,1])
plot(linspace(0,1,100000), yy, 'LineWidth', 2.0);

subplot(2,3,2)
xlim([0,1])
hold on
ylim([0,1])
plot(linspace(0,1,100000), yy, 'LineWidth', 2.0);

subplot(2,3,3)
xlim([0,1])
hold on
ylim([0,1])
plot(linspace(0,1,100000), yy, 'LineWidth', 2.0);

%% sobol sequence
for idx = 1:numMax
    subplot(2,3,1)
    scatter(ss(idx), 0, 'black', 'filled');
    hold on
    subplot(2,3,4)
    hold on
    errr = abs(sum(funct(ss(1:idx)))/idx - truth);
    scatter(idx, errr, 'red', 'filled')
    hold on
    xlim([0,numMax])
end

%% uniform sampling
for idx = 1:numMax
    subplot(2,3,2)
    unis(idx) = rand(1);
    scatter(unis(idx), 0, 'black', 'filled');
    hold on
    subplot(2,3,5)
    hold on
    errr = abs(sum(funct(unis(1:idx)))/idx - truth);
    scatter(idx, errr, 'red', 'filled')
    hold on
    xlim([0,numMax])
end
%% randomized QMC sampling
ss = scramble(ss,'MatousekAffineOwen');
for idx = 1:numMax
    subplot(2,3,3)
    scatter(ss(idx), 0, 'black', 'filled');
    hold on
    subplot(2,3,6)
    hold on
    errr = abs(sum(funct(ss(1:idx)))/idx - truth);
    scatter(idx, errr, 'red', 'filled')
    hold on
    xlim([0,numMax])
end

function val = funct(in)
    val = normpdf(in, 0.4, 0.005)/100;
end