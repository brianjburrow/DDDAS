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

for iCenter = 1:numSamples
    currentMean = means(:, iCenter);
    %disp(currentMean)
    sampleSet(iCenter).samples = mvnrnd(currentMean, Sig_max, numInnerSamples)';
    %sampleSet(iCenter).samples = mvnrnd_truncated(currentMean, Sig_max, numInnerSamples);
    scatter(sampleSet(iCenter).samples(1,:), sampleSet(iCenter).samples(2,:), 5, 'filled', 'blue');
    hold on
end

hold on
scatter(means(1,:), means(2,:), 10, 'filled', 'black')
hold on
% xlim([0,1])
% ylim([0,1])



function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'red', 'LineWidth', 2.0);
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