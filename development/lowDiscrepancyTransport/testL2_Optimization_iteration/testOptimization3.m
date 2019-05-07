%% Generate GMD
clear all
close all
clc

numSamples = 2000;
numInnerSamples = 100;

%means = normrnd(0, 1, [2, numSamples]).^2;
means = random_sampler_joint(numSamples)';
%means = damage_example(numSamples);

%Sig_max = 0.05*[1, 0; 0, 1];
Sig_max = diag(var(means'));


N = 1000;

yy = @(xx) (log(N).^xx)/N;
zz = @(xx) 1./sqrt(xx);


figure(1)
subplot(2,2,1)
for iCenter = 1:numSamples
    currentMean = means(:, iCenter);
    if iCenter == 1
        samples = mvnrnd(currentMean, Sig_max, numInnerSamples)';
    else
        samples = [samples, mvnrnd(currentMean, Sig_max, numInnerSamples)'];
    end
end
scatter(samples(1,:), samples(2,:), 10, 'filled')
hold on
scatter(means(1,:), means(2,:), 10, 'filled')
xlim([-1, 10])
ylim([-1, 10])

covSet = zeros(length(means'), 2,2);

sigMax = cov(samples');
for iCount = 1:length(means')
    covSet(iCount,:,:) = sigMax/20;
end

tic
log_pSet = multiVariateMixtureLogPdf(samples', means', covSet);
disp(toc)

bounds = [min(samples'); max(samples')];

largeFunc = @(xx) multiVariateMixtureLogPdf(xx, means', covSet);
smallFunc = @(xx) multiVariateMixtureLogPdf(xx, means', covSet);

sobolOptions.maxCardinality = 10000;
sobolOptions.thresholds = 0.05;

tic
[library, indexSet] = constructLibrary(largeFunc, ...
    smallFunc, ...
    sobolOptions,...
    bounds);
disp(toc)
xlim([min(library(:,1)), max(library(:,1))])
ylim([min(library(:,2)), max(library(:,2))])
subplot(2,2,4)
scatter(library(:,1), library(:,2), 10, 'filled')
xlim([min(library(:,1)), max(library(:,1))])
ylim([min(library(:,2)), max(library(:,2))])

figure(2)
ss = sobolset(2);
tsamples = [linspace(-1,10,1000)', linspace(-1,10,1000)'];
xi = convertSobolToBounds(tsamples, bounds);
f = largeFunc(xi);
[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);
hold on
scatter(xi(:,1), xi(:,2), 2.5, 'filled')
hold on
[C,h] = contour(X,Y,Z, 20);
colormap(hot)
xlim([min(library(:,1)), max(library(:,1))])
ylim([min(library(:,2)), max(library(:,2))])

function p = multiVariateMixtureLogPdf(xx, muSet, covSet)
    % xx is a nSample x nDim matrix
    % muSet is an nDist x nDim matrix
    % covSet is an nDim x nDim x nDist Array
    % returns p: an nSample x 1 array
    
    [nDist, ~] = size(muSet);
    [nSamples, nDim] = size(xx);
    p = zeros([nSamples, nDist]);
    for iDist = 1:nDist
        p(:, iDist) = mvnpdf(xx, muSet(iDist, :), reshape(covSet(iDist, :,:), nDim, nDim));
    end
    p = max(p, [], 2);
end

function samples = convertSobolToBounds(samples, bounds)
        % bounds =  2 x nDim matrix, lower bound is first row, upper = 2nd
        [nSamp, nDim] = size(samples);
        blankOnes = ones([nSamp, 1]);
        for iDim = 1:nDim
            samples(:,iDim) = samples(:,iDim)*(...
                bounds(2,iDim) - bounds(1,iDim)...
                ) + ...
                bounds(1,iDim)*blankOnes;
        end
end
    
function theta = random_sampler_joint(numSamples)
    %theta = zeros([2,numSamples]);
    r_1 = normrnd(0, 1, [1, numSamples]);
    r_2 = normrnd(0, 1, [1, numSamples]);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end

function SAMPLES = damage_example(numSamples)
    nTime = 500;
    nTrials = numSamples/nTime;
    samples(nTrials).samps = zeros(nTime, 2);
    parfor iTrial = 1:nTrials 
        damage = [0,0];
        for iTime = 1:nTime
            add = mvnrnd_truncated([0,0], 0.01*eye(2), 1)';
            samples(iTrial).samps(iTime, :) = damage + add;
        end
    end
    SAMPLES = samples(1).samps;
    for iTrial = 2:nTrials
        SAMPLES = [SAMPLES; samples(iTrial).samps];
    end
    
    SAMPLES = SAMPLES';
end

function samples = mvnrnd_truncated(mean, coVar, numSamples)
    % Bound to 2-D hypercube
    for iSamp = 1:numSamples
        continueLoop = 1;
        counter = 0;
        while continueLoop
            prop = mvnrnd(mean, coVar);
            if (prop(1,1) < 0) || (prop(1,1) > 1) || (prop(1,2) < 0) || (prop(1,2) > 1)
                counter = counter + 1;
                continue
            else
                samples(iSamp,:) = prop;
                continueLoop = 0;
            end
        end
    end
    samples = samples';
end