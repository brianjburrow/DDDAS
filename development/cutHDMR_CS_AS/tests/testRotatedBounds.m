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
cutLine = [-.25, .75];                                                           % The cutLine about which cutHDMR is computed (reference point essentially)

%% Create Active Subspace
funcHandle = @(xx) testFunc(xx);
gradHandle = @(xx) gradTestFunc(xx);

actSub = ActiveSubspace(gradHandle, ...
    nMCsamples, ...
    nDimensions,...
    samplerHandle, ...                                                      % BB May 2, 2019: should this be a different random variable than chebyshev?
    eigCut);                                                                % Cumulative eigenvalues is a good way to cut off dimensions


%% Create altered cutHDMR surrogate by projecting onto active subspace, computation accelerated using Compressed Sensing
% Map all necessary components to the active subspace
funcHandle2 = @(yy) testFuncMapped(yy, actSub.W1);                          % Reparameterize the test function to take values on the mapped subspace.
samplerHandle2 = @(numSamp) actSub.ProjectSamples(samplerHandle(numSamp));  % Project chebyshev points onto the new subspace.
cutLine2 = actSub.ProjectSamples(cutLine);      

subplot(2,2,1)
samples1 = samplerHandle(1000);
scatter(samples1(:,1), samples1(:,2), 10, 'k', 'filled')
hold on
xlim([-2, 2])
ylim([-2, 2])
xlabel('X_1')
ylabel('X_2')


subplot(2,2,2)
samples2 = actSub.ProjectSamples(samples1);
[sampCut1, sampCut2, orig1, orig2] = rejectionSamplerActSub(1000, samplerHandle, actSub, [-1, -1], [1, 1], cutLine);
scatter(samples2(:,1), samples2(:,2), 10, 'b', 'filled')
hold on
scatter(orig1(:,1), orig1(:,2), 10, 'r', 'filled')
scatter(orig2(:,1), orig2(:,2), 10, 'r', 'filled')
xlim([-2, 2])
ylim([-2, 2])
xlabel('Y_1')
ylabel('Y_2')

lowerBounds = [-1; -1];
upperBounds = [1; 1];
rotationMatrix = actSub.W1';
[newLowBounds, newUpBounds] = FindRotatedDomain(...
    lowerBounds, upperBounds, rotationMatrix...
    );

subplot(2,2,3)
scatter(samples1(:,1), samples1(:,2), 10, 'k', 'filled')
hold on
scatter(sampCut1(:,1), sampCut1(:,2), 10, 'r', 'filled')
scatter(sampCut2(:,1), sampCut2(:,2), 10, 'r', 'filled')
xlim([-2, 2])
ylim([-2, 2])
xlabel('X_1')
ylabel('X_2')


subplot(2,2,4)
z = ConvertRotatedDomain(samples2, newLowBounds, newUpBounds);
scatter(z(:,1), z(:,2), 10, 'b', 'filled');
hold on
xlim([-2, 2])
ylim([-2, 2])
xlabel('Y_1^{u}')
ylabel('Y_2^{u}')

orig1 = ConvertRotatedDomain(orig1, newLowBounds, newUpBounds);
orig2 = ConvertRotatedDomain(orig2, newLowBounds, newUpBounds);
scatter(orig1(:,1), orig1(:,2), 10, 'r', 'filled')
scatter(orig2(:,1), orig2(:,2), 10, 'r', 'filled')

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

function [samp1, samp2, orig1, orig2] = rejectionSamplerActSub(nSamp, sampHand, acSubspace, lbs, ubs, cutLine)
    orig1 = [];
    orig2 = [];
    samp1 = [];
    samp2 = [];
    loop = 1;
    while loop
        samples1 = sampHand(nSamp);
        samples2 = acSubspace.ProjectSamples(samples1);
        samp3 = acSubspace.invertProjection([samples2(:,1), cutLine(2)*ones(size(samples2(:,1)))]);     % Get cut samples in new subspace
        samp4 = acSubspace.invertProjection([cutLine(1)*ones(size(samples2(:,2))), samples2(:,2)]);     % get cut samples in new subspace
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
            samp3(list, :) = [];                                            % Remove rows
            samp1 = [samp1; samp3];
            samp4(list, :) = [];
            samp2 = [samp2; samp4];
            nSamp = nSamp - length(samp2(:,1));
            samp = [samples2(:,1), cutLine(2)*ones(size(samples2(:,1)))];
            samp(list, :) = [];
            orig1 = [orig1; samp];
            samp = [cutLine(1)*ones(size(samples2(:,2))), samples2(:,2)];
            samp(list, :) = [];
            orig2 = [orig2; samp];
        end
    end
end
