clear all
close all
clc
%% testDynamicActiveSubspace.m
% Author: Brian Burrows
% Computes a series of Active Subspaces over the same function given a
% series of dynamically evolving probability distributions (i.e., 
% what you would encounter in filtering).

%% Test sampling distribution 1
ss = sobolset(2);
%randomSampler = @(nSamp) ss(1:nSamp, :)*15 + [-5*ones([nSamp,1]), zeros([nSamp, 1])];
nDimensions = 2;
nSamples = 100;                                                                  % Number of training samples for the surrogate
nMCsamples = 10000;                                                              % Number of Monte carlo samples used to train the active subspace


%% Create Active Subspaces
figure(1)
title("Dynamically Evolving Active Subspace")
hold on
xlim([-5, 10])
ylim([0, 15])
nSamp = 40;
muSet = ss(1:nSamp, :)*15 + [-5*ones([nSamp,1]), zeros([nSamp, 1])];
cov = 0.25*eye(2);

[nMu, ~] = size(muSet);
for iMu = 1:nMu
    randomSampler = @(nSamples) mvnrnd(muSet(iMu,:), cov, nSamples);
    eigCut = 0.985;
    braninHandle = @(xx) BraninGradient(xx);
    samplerHandle = @(nSamp) randomSampler(nSamp);
    actSub = ActiveSubspace(braninHandle, ...
        nMCsamples, ...
        nDimensions,...
        samplerHandle, ...
        eigCut);  % Cumulative eigenvalues is a good way to cut off dimensions
    l1 = quiver(muSet(iMu, 1), muSet(iMu, 2), actSub.W1(1,1), actSub.W1(2,1), 'k' ,'DisplayName','Local Basis', 'MaxHeadSize', 1);
    try
        quiver(muSet(iMu, 1), muSet(iMu, 2), actSub.W1(1,2), actSub.W1(2,2), 'k','DisplayName','Local Basis', 'MaxHeadSize', 1)
    catch
        fprintf("Second active basis vector does not exist in this location \n")
    end
end

legend([l1])
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