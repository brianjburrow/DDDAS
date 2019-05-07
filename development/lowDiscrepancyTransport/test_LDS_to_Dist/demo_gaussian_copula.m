%% demo_gaussian_copula
clear all
close all
clc

%% Generate samples from an uncorrelated Gaussian random variable

n = 500;                                                                    % select number of samples
Z = mvnrnd([0 0], [1 0; 0 1], n);                                           % generate using MATLAB's multivariate random normal sampler

figure(1)
scatter(Z(:,1), Z(:,2))
hold on
title("MVNRND")
%% Generate samples from an uncorrelated Gaussian random variable
u = rand([n, 2]);                                                           % generate uniform samples using MATLAB's "rand" function
X = uniformTOgaussian(u);                                                   % use inverse transform method to generate standard gaussian samples

figure(2)
scatter(X(:,1), X(:,2))
hold on
title("rand: INV Transform")
%% Generate low-discrepancy points from uncorrelated Gaussian random variable
ss = sobolset(2);                                                           % initialize sobol set in two dimensions
u = ss(1:n, :);                                                             % generate deterministic samples from sobol set that appear uniformly distributed
Y = uniformTOgaussian(u);                                                   % use inverse transform to generate deterministic samples that appear standard Gaussian
figure(3)
scatter(Y(:,1), Y(:,2))
hold on
title("sobol: INV Transform")

figure(2)
scatter(X(:,1), X(:,2))
hold on
title("rand: INV Transform")
%% Generate low-discrepancy points from uncorrelated Gaussian random variable
ss = sobolset(2);                                                           % initialize sobol set in two dimensions
ss = scramble(ss, 'MatousekAffineOwen');
u = ss(1:n, :);                                                             % generate deterministic samples from sobol set that appear uniformly distributed
Y = uniformTOgaussian(u);                                                   % use inverse transform to generate deterministic samples that appear standard Gaussian
figure(4)
scatter(Y(:,1), Y(:,2))
hold on
title("randomized sobol: INV Transform")

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end