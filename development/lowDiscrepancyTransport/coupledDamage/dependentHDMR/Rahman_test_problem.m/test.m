%% This file is meant to get some visual insight into an example from
% "A generalized ANOVA dimensional decomposition for dependent probability
% measures." by Sharif Rahman

%% Case 1, uncorrelated Gaussian random variables
close all

mu1 = [0, 0, 0];
cov1 = [1, 0, 0;
        0, 1, 0;
        0, 0, 1];
figure(1)
plotPairwiseMarginals(mu1, cov1)
suptitle("Uncorrelated Inputs")
%% Case 2, equally correlated Gaussian random variables

mu2 = [0, 0, 0];
cov2 = [1, 1/5, 1/5;
        1/5, 1, 1/5;
        1/5, 1/5, 1];
figure(2)
plotPairwiseMarginals(mu2, cov2)
suptitle("Equally Correlated Inputs")

%% Case 3, unequally correlated Gaussian random variables

mu3 = [0, 0, 0];
cov3 = [1, 1/5, 2/5;
        1/5, 1, 4/5;
        2/5, 4/5, 1];
figure(3)
plotPairwiseMarginals(mu3, cov3)
suptitle("Unequal, Positively Correlated Inputs")


%% Case 4, unequally correlated Gaussian random variables

mu4 = [0, 0, 0];
cov4 = [1, -1/5, 2/5;
        -1/5, 1, -4/5;
        2/5, -4/5, 1];
figure(4)
plotPairwiseMarginals(mu4, cov4)
suptitle("Unequal, Mixed Correlated Inputs")

%% Define function of interest

function y = testFunc(X)
    x1 = X(:, 1);
    x2 = X(:, 2);
    x3 = X(:, 3);
    y = 12*ones(size(x1)) + 4*x1 + 4*x2 + 4*x3 + ...
        x1.*x2 + x1.*x3 + x2.*x3;
end


function plotPairwiseMarginals(muVec, covMat)
    counter = 1;
    for iFeature = 1:2
        for iContrast = iFeature+1:3
            subplot(2, 2, counter)
            hold on
            tempMu = [muVec(iFeature), muVec(iContrast)];
            tempCov = [covMat(iFeature, iFeature), covMat(iFeature, iContrast);
                covMat(iFeature, iContrast), covMat(iContrast, iContrast)];
            counter = counter + counter;
            samples = mvnrnd(tempMu, tempCov, 10000);
            scatter(samples(:,1), samples(:,2), 10, 'filled', 'k');
            xlabel(sprintf('%d', iFeature))
            ylabel(sprintf('%d', iContrast))
        end
    end
end