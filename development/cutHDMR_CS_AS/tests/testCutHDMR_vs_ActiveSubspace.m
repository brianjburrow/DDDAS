clear all
close all
clc
%% Setup problem 
nMCsamples = 10000;                                                         % Number of Monte Carlo samples used to construct active subspace
nDimensions = 2;                                                            % Number of input variables into our function of interest
eigCut = 0.99999999;                                                        % Truncate dimension when the cumulative sum of the eigenvalues are greater than eigCut
samplerHandle = @(numSamp) [...                                             % A random variable generator over the input space for constructing the HDMR and Active subspace.
    chebyshevSampler(numSamp, -1, 1),...
    chebyshevSampler(numSamp, -1, 1)...
    ];
cutLine = [0, 0];                                                           % The cutLine about which cutHDMR is computed (reference point essentially)

%% Create Active Subspace
funcHandle = @(xx) testFunc(xx);
gradHandle = @(xx) gradTestFunc(xx);

actSub = ActiveSubspace(gradHandle, ...
    nMCsamples, ...
    nDimensions,...
    samplerHandle, ...                                                      % BB May 2, 2019: should this be a different random variable than chebyshev?
    eigCut);                                                                % Cumulative eigenvalues is a good way to cut off dimensions


%% Create Unaltered cutHDMR surrogate, computation accelerated using Compressed Sensing and convergence monitoring
nSamples = 3;                                                               % Number of initial Samples to evaluate 
hdmr = cutHDMR_compressedSensing_withActiveSubspaces(funcHandle, cutLine, nSamples, 4,...       % Initialize cutHDMR object
    'legendre', samplerHandle);
hdmr = hdmr.run();                                                          % Construct cutHDMR surrogate model

%% Create altered cutHDMR surrogate by projecting onto active subspace, computation accelerated using Compressed Sensing
% Map all necessary components to the active subspace
funcHandle2 = @(yy) testFuncMapped(yy, actSub.W1);                          % Reparameterize the test function to take values on the mapped subspace.
samplerHandle2 = @(numSamp) actSub.ProjectSamples(samplerHandle(numSamp));  % Project chebyshev points onto the new subspace.
cutLine2 = actSub.ProjectSamples(cutLine);                                  % Project the cutpoints onto the new subspace.

figure(100)
samp1 = samplerHandle(1000);
samp2 = actSub.ProjectSamples(samp1);
scatter(samp1(:,1), samp1(:,2), 10, 'k', 'filled')
hold on
scatter(samp2(:,1), samp2(:,2), 10, 'b', 'filled')
figure(1)

cutSamples = [zeros(size(samp2(:,1))), samp1(:,2)];
cutSamp2 = actSub.ProjectSamples(cutSamples);
scatter(cutSamp2(:,1), cutSamp2(:,2))

figure(10000)
whatSamples = actSub.invertProjection(samp1);
scatter(whatSamples(:,1), whatSamples(:,2))

aoeue

nSamples = 3;                                                               % Number of initial Samples to evaluate for cutHDMR
hdmr2 = cutHDMR_compressedSensing(funcHandle2, cutLine2, nSamples, 4,...
    'legendre', samplerHandle2);
hdmr2 = hdmr2.run();

samp1 = [linspace(-1,1)', linspace(-1,1)'];
samp2 = actSub.ProjectSamples(samp1);

[XX1, YY1] = meshgrid(samp1(:,1)', samp1(:,2)');

[XX2, YY2] = meshgrid(samp2(:,1)', samp2(:,2)');

for iY = 1:100
    for iX = 1:100
        ZZ1(iX, iY) = testFunc([XX1(iX, iY), YY1(iX, iY)]);
        ZZ2(iX, iY) = testFuncMapped2([XX1(iX, iY), YY1(iX, iY)], actSub.W1);
    end
end

n = get(gcf, 'Number');
figure(n + 1)
surf(XX1, YY1, ZZ1)
figure(n + 2)
surf(XX1, YY1, ZZ2)

function output = testFunc(X)
    output = 1 + X(:,1) + 0.5*(3*X(:,2).^2 - 1) + X(:,1).*X(:,2);
end

function output = testFuncMapped(Y, W)
    X = (W * Y')';
    output = 1 + X(:,1) + 0.5*(3*X(:,2).^2 - 1) + X(:,1).*X(:,2);
end

function output = testFuncMapped2(X, W)
    Y = (W' * X')';
    X = (W * Y')';
    output = 1 + X(:,1) + 0.5*(3*X(:,2).^2 - 1) + X(:,1).*X(:,2);
end

function gradMatrix = gradTestFunc(X)
    [nSamples, nDim] = size(X);
    gradMatrix = zeros([nDim, nSamples]);
    gradMatrix(1,:) = ones([1, nSamples]) + X(:,2)';
    gradMatrix(2,:) = 3*X(:,2)' + X(:,1)';
end

function [samp1, samp2] = rejectionSamplerActSub(nSamp, sampHand, acSubspace, lbs, ubs)
    samp1 = [];
    samp2 = [];
    loop = 1;
    while loop
        samples1 = sampHand(nSamp);
        samples2 = acSubspace.ProjectSamples(samples1);
        samp3 = acSubspace.invertProjection([samples2(:,1), zeros(size(samples2(:,1)))]);     % Get cut samples in new subspace
        samp4 = acSubspace.invertProjection([zeros(size(samples2(:,2))), samples2(:,2)]);     % get cut samples in new subspace
        nDim = length(lbs);
        list = [];
        for iSamp = 1:nSamp
            loop = 1;
            for iDim = 1:nDim
                if (samp3(iSamp, iDim) < lbs(iDim)) || (samp3(iSamp, iDim) > ubs(iDim))
                    list = [list; iSamp];
                elseif (samp4(iSamp, iDim) < lbs(iDim)) || (samp4(iSamp, iDim) > ubs(iDim))
                    list = [list; iSamp];
                end
            end
        end
        if isempty(list)
            loop = 0;
        else
            samp3(list, :) = [];                                           % Remove rows
            samp1 = [samp1; samp3];
            samp4(list, :) = [];
            samp2 = [samp2; samp4];
            nSamp = nSamp - length(samp2(:,1));
        end
    end
end