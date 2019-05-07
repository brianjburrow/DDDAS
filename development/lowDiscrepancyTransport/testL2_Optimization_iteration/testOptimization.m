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


% sampWhat = mvnrnd_truncated([0.002, 0.25], Sig_max, numInnerSamples)';
% scatter(sampWhat(:,1), sampWhat(:,2), 5, 'black')
% hold on
% scatter(0.002, 0.25, 10, 'red')
% hold on

N = 1000;

yy = @(xx) (log(N).^xx)/N;
zz = @(xx) 1./sqrt(xx);



for iCenter = 1:numSamples
    currentMean = means(:, iCenter);
    %disp(currentMean)
    sampleSet(iCenter).samples = mvnrnd(currentMean, Sig_max, numInnerSamples)';
    %sampleSet(iCenter).samples = mvnrnd_truncated(currentMean, Sig_max, numInnerSamples);
    scatter(sampleSet(iCenter).samples(1,:), sampleSet(iCenter).samples(2,:), 5, 'filled', 'blue');
    hold on
end

scatter(means(1,:), means(2,:), 10, 'filled', 'black')
xlim([0,1])
ylim([0,1])
numSob = 5000;
ss = sobolset(2);
pointSet = ss(1:numSob, :);
for idx = 1:numSob
    scatter(pointSet(idx,1), pointSet(idx,2), 5*rand(1), 'filled', 'red')
    hold on
end


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

function samples = damage_example(numSamples)
    nTime = 100;
    nTrials = numSamples/nTime;
    counter = 1;
    for iTrial = 1:nTrials 
        damage = [0.25,0.4];
        for iTime = 1:nTime
            add = mvnrnd_truncated([0,0], 0.005*eye(2), 1)';
            samples(counter, :) = damage + add;
            counter = counter + 1;
        end
    end
    samples = samples';
end

function samples = mvnrnd_truncated(mean, coVar, numSamples)
    % Bound to 2-D hypercube
    for iSamp = 1:numSamples
        continueLoop = 1;
        counter = 0;
        while continueLoop
            prop = mvnrnd(mean, coVar);
            if (prop(1,1) < 0.25) || (prop(1,1) > 1) || (prop(1,2) < 0.4) || (prop(1,2) > 1)
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