clear all
close all
clc

%% testGaussianActiveSubspace.m
% Author: Brian Burrows
% Constructs surrogates over the full parameter space using quasi-random
% uniform samples to compare the accuracy loss of the surrogate when
% going from the full problem dimension (2), to a lower dimension (depends
% on the value for eigCut).  I compare two approaches, a surrogate trained on
% the active subspace (either 1 or 2 D depending on eigCut) and a surrogate
% trained on samples from the original parameter space.


%% Test sampling distribution 1
ss = sobolset(2);
%randomSampler = @(nSamples) mvnrnd([3,7], 4*eye(2), nSamples);
randomSampler = @(nSamp) ss(1:nSamp, :)*15 + [-5*ones([nSamp,1]), zeros([nSamp, 1])];
nDimensions = 2;
nSamples = 100;                                                                  % Number of training samples for the surrogate
nMCsamples = 500;                                                                % Number of Monte carlo samples used to train the active subspace


%% test the functions
testRVS = randomSampler(nSamples);                                               % Test the random sampler function to make sure it returns the correct shape
testBranin = BraninGradient(testRVS);
quiver(testRVS(:,1)', testRVS(:,2)', testBranin(1,:), testBranin(2,:))
hold on
xlim([-5, 10]);
ylim([0, 15]);
%% Create Active Subspace
eigCut = 0.99;
braninHandle = @(xx) BraninGradient(xx);
samplerHandle = @(nSamp) randomSampler(nSamp);
actSub = ActiveSubspace(braninHandle, ...
    nMCsamples, ...
    nDimensions,...
    samplerHandle, ...
    eigCut);  % Cumulative eigenvalues is a good way to cut off dimensions

%% Project Input Samples onto Active Subspace
projTestRVS = actSub.ProjectSamples(testRVS);
outputs = BraninFunction(testRVS);                                          % Must use full samples

%% Fit surrogate
figure(2)
gp = fitrgp(projTestRVS, outputs);

x = linspace(-5, 10)';
y = linspace(0, 15)';

ZZZ1 = zeros(100);
XXX = zeros(100);
YYY = zeros(100);
for iX = 1:100
    for iY = 1:100
        XXX(iX, iY) = x(iX);
        YYY(iX, iY) = y(iY);
        ZZZ1(iX, iY) = predict(gp, actSub.ProjectSamples([x(iX), y(iY)]));
    end
end

surf(XXX, YYY, ZZZ1);
hold on
scatter3(testRVS(:,1), testRVS(:,2), outputs, 'k', 'filled') 
title("Active Subspace Surrogate")

figure(3)

gp2 = fitrgp(testRVS, outputs);
ZZZ2 = zeros(100);
for iX = 1:100
    for iY = 1:100
        ZZZ2(iX, iY) = predict(gp2, [x(iX), y(iY)]);
    end
end
surf(XXX, YYY, ZZZ2);
hold on
scatter3(testRVS(:,1), testRVS(:,2), outputs, 'k', 'filled') 
title("Normal Subspace Surrogate")

figure(4)
ZZZ3 = zeros(100);
for iX = 1:100
    for iY = 1:100
        ZZZ3(iX, iY) = BraninFunction([x(iX), y(iY)]);
    end
end
surf(XXX, YYY, ZZZ3);
hold on
title("Actual Function")

figure(5)
surf(XXX, YYY, abs(ZZZ1 - ZZZ3))
hold on
fprintf('MSE Active Subspace: %f \n', mean(mean((ZZZ1 - ZZZ3).^2)))
title("Active Subspace Surrogate Error")

figure(6)
surf(XXX, YYY, abs(ZZZ2 - ZZZ3))
fprintf('MSE Normal Subspace: %f \n', mean(mean((ZZZ2 - ZZZ3).^2)))
title("Normal Subspace Surrogate Error")

%% See what mapped distribution looks like after transforming to active subspace
n = get(gcf, 'Number')
figure(n + 1)
samps = randomSampler(1000);
newSamps = actSub.ProjectSamples(samps);
subplot(2,2,1)
scatter(samps(:,1), samps(:,2), 10, 'k', 'filled')
subplot(2,2,4)
scatter(newSamps(:,1), newSamps(:,2), 10, 'k', 'filled')

function output = BraninFunction(X)
    % https://www.sfu.ca/~ssurjano/branin.html
    x1 = X(:,1);
    x2 = X(:,2);
    a = 1;
    b = 5.1 / (4*pi^2);
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    [nSamples, ~] = size(x1);
    oneVec = ones(nSamples, 1);
    output = a * (x2 - b * x1.^2 + c*x1 - r*oneVec).^2 + ...
        s*(1 - t)*cos(x1) + s*oneVec;
end

function gradMatrix = BraninGradient(X)
    x1 = X(:,1);
    x2 = X(:,2);
    a = 1;
    b = 5.1 / (4*pi^2);
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    [nSamples, ~] = size(x1);
    oneVec = ones(nSamples, 1);
    % partial calculation
    repeatedCalc = (x2 - b*x1.^2 + c*x1 - r*oneVec);
    % component 1 
    part1 = (2*a*repeatedCalc.*(-2*b*x1 + c*oneVec) - s*(1 - t)*sin(x1))';
    part2 = (2*a*repeatedCalc)';
    gradMatrix = [part1 ; part2];
end