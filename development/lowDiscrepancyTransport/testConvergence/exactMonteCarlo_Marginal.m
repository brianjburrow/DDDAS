clear all
clc
close all
rng(0)
numRepeats = 100;
sampleSize = 25:25:400;

%% Monte Carlo Analysis
numTargets = 500;
figure(100)
xx = random_sampler(1000000);
subplot(2,1,1)
scatter(xx(:,1), xx(:,2))
yy = forwardModel(xx);
subplot(2,1,2)
scatter(xx(:,1), yy)
truth = mean(forwardModel(random_sampler(100000)));
%for iRepeat = 1:numRepeats
iRepeat = 1;
samples = random_sampler(numTargets);
    for iSampSize = 1:length(sampleSize)
        outputs = forwardModel(samples);
        moment1(iRepeat, iSampSize) = abs(mean(outputs) - truth);
    end
%end

figure(1)
subplot(1,1,1)
plot(sampleSize, moment1, 'black', 'LineWidth', 2.0, 'DisplayName', 'Target Samples')
hold on
xlabel("Sample Size")
ylabel("Absolute Error")
drawnow
moment12 = zeros(numRepeats, length(sampleSize));                           % uniform random: numRepeats
moment13 = zeros(numRepeats, length(sampleSize));                           % joint random:   numRepeats
moment14 = zeros(1, length(sampleSize));                                    % Sobol uniform:  deterministic so no repeats
moment15 = zeros(1, length(sampleSize));                                    % Sobol joint:    deterministic so no repeats
%% Uniform Random Sampling with L2O Analysis
for iRepeat = 1:numRepeats
    disp(iRepeat)
    tic
    parfor iSampSize = 1:length(sampleSize)
        %samples = random_sampler(sampleSize(iSampSize));
        % Proposal Samples
        propSamples = 10*rand(sampleSize(iSampSize), 2);
%         
%         [indx, ~] = knnsearch(propSamples, samples,'K', 1);
%         indx = sort(indx);
%         propSamples = propSamples(indx,:); 
         [weights, propSamples, EFF] = determineWeights(propSamples, samples);
        %disp(EFF)
        outputs = forwardModel(propSamples);

        moment12(iRepeat,iSampSize) = abs(weightedMean(outputs, weights) - truth);
    end
    toc
end

figure(1)
subplot(1,1,1)
plot(sampleSize, mean(moment12), 'red', 'LineWidth', 2.0, 'DisplayName', 'Uniform Samples')
hold on
%errorbar(sampleSize, mean(moment12), 3*std(moment12), 'blue')
hold on
xlabel("Sample Size")
ylabel("Absolute Error")
drawnow

figure(2)
subplot(1,1,1)
plot(sampleSize, var(moment12), 'red', 'LineWidth', 2.0, 'DisplayName', 'Uniform Samples')
hold on
%errorbar(sampleSize, mean(moment12), 3*std(moment12), 'blue')
hold on
xlabel("Sample Size")
ylabel("Estimate Variance")
drawnow

disp('a')

%% Joint Random Sampling with L2O Analysis
for  iRepeat = 1:numRepeats
    disp(iRepeat)
    parfor iSampSize = 1:length(sampleSize)
        %samples = random_sampler(sampleSize(iSampSize));
        % Proposal Samples
        propSamples = random_sampler_joint(sampleSize(iSampSize));
%         [indx, ~] = knnsearch(propSamples, samples,'K', 1);
%         indx = sort(indx);
%         propSamples = propSamples(indx,:); 
% 
         [weights, propSamples, EFF] = determineWeights(propSamples, samples);
%         disp(EFF)

        outputs = forwardModel(propSamples);
        moment13(iRepeat, iSampSize) = abs(weightedMean(outputs, weights)- truth);
    end
end

figure(1)
subplot(1,1,1)
plot(sampleSize, mean(moment13), 'blue', 'LineWidth', 2.0, 'DisplayName', 'Joint Samples')
hold on
%errorbar(sampleSize, mean(moment13), 3*std(moment12), 'blue')
hold on
legend()
xlabel("Sample Size")
ylabel("Absolute Error")
drawnow

figure(2)
subplot(1,1,1)
plot(sampleSize, var(moment13), 'blue', 'LineWidth', 2.0, 'DisplayName', 'Uniform Samples')
hold on
%errorbar(sampleSize, mean(moment12), 3*std(moment12), 'blue')
hold on
xlabel("Sample Size")
ylabel("Estimate Variance")

aoeueo
%% Uniform Low Discrepancy Sequence with L2O Analysis
for iRepeats = 1:numRepeats
    parfor iSampSize = 1:length(sampleSize)
        %samples = random_sampler(sampleSize(iSampSize));
        % Proposal Samples
        ss = sobolset(2)
        ss = scramble(ss, 'MatousekAffineOwen');
        propSamples = 10*ss(1:sampleSize(iSampSize), :);
%         [indx, ~] = knnsearch(propSamples, samples,'K', 1);
%         indx = sort(indx);
%         propSamples = propSamples(indx,:); 
        
        [weights, propSamples, EFF] = determineWeights(propSamples, samples);
    %    disp(EFF)
        outputs = forwardModel(propSamples);
        moment14(1,iSampSize) = abs(weightedMean(outputs, weights) - truth);
    end
end
figure(1)
plot(sampleSize, moment14, 'green', 'LineWidth', 2.0, 'DisplayName', 'Sobol Samples')
hold on
%% Joint Low Discrepancy Sequence with L2O Analysis

for iRepeats = 1:numRepeats
    parfor iSampSize = 1:length(sampleSize)
        %samples = random_sampler(sampleSize(iSampSize));
        % Proposal Samples
        propSamples = random_sampler_2(sampleSize(iSampSize));
        [indx, ~] = knnsearch(propSamples, samples,'K', 1);
        indx = sort(indx);
        propSamples = propSamples(indx,:); 
        
        [weights, propSamples, EFF] = determineWeights(propSamples, samples);
    %         disp(EFF)
        outputs = forwardModel(propSamples);
        moment15(1, iSampSize) = abs(weightedMean(outputs, weights)- truth);
    end
end

figure(1)
subplot(1,1,1)
plot(sampleSize, moment15, 'cyan', 'LineWidth', 2.0, 'DisplayName', 'Sobol Joint Samples')
hold on
legend()
xlabel("Sample Size")
ylabel("Absolute Error")
drawnow

%% Begin function definitions
function theta = random_sampler(numSamples)
    %theta = zeros([2,numSamples]);
    r_1 = normrnd(1, 0.5, [1, numSamples]);
    r_2 = normrnd(1, 0.5, [1, numSamples]);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end

function output = forwardModel(x)
    output = x(:,1).^2 + cos(x(:,2));
end

function theta = random_sampler_joint(numSamples)
    %theta = zeros([2,numSamples]);
    r_1 = normrnd(0, 1, [1, numSamples]);
    r_2 = normrnd(0, 1, [1, numSamples]);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end

function [weights, samples, EFF] = determineWeights(samples, newSamples)
    %% Determine Weights
    % samples are proposal samples
    % newSamples are target samples
    itr = 20000;
    samples = samples';
    newSamples = newSamples';
    numSamples = length(samples(1,:));
    concatSamps = [samples, newSamples];
    mins = min([samples, newSamples], [], 2);
    maxs = max([samples, newSamples], [], 2);
    normedProps = (samples - mins)./(maxs - mins);
    normedTargets = (newSamples - mins)./(maxs - mins);
    [normedProps, propIdx] = sortrows(normedProps');
    [normedTargets, ~] = sortrows(normedTargets');
    samples = samples(:,propIdx);

    w = ones(numSamples,1)./numSamples;
    H = matvecProd(normedProps, w);
    B = buildXYb(normedProps, normedTargets);
    % Construct B vector
    feval = zeros(itr+1, 1);
    [weights, feval] = frankwolfe(normedProps, B, w, H, feval, itr);
     EFF = sum(weights)^2 / (sum(weights.^2));
     samples = samples';
end
function theta = random_sampler_2(numSamples)
    %theta = zeros([2,numSamples]);
    ss = sobolset(2);
    ss = scramble(ss, 'MatousekAffineOwen');
    uniSamples = ss(1:numSamples, :)';
    refSamples = uniformTOgaussian(uniSamples);
    r_1 = refSamples(1,:);
    r_2 = refSamples(2,:);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end

function meanVal = weightedMean(weights, outputs)
    numSamp = length(outputs(:,1));
    numDim = length(outputs(1,:));
    meanVal = zeros(1, numDim);
    for iSamp = 1:length(weights)
        meanVal(1,:) = meanVal(1,:) + weights(iSamp) * outputs(iSamp,:);
    end
end

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end
