clear all
close all
clc
%% Set up test parameters
nTrainSamples = 15;                                                         % Original size of design of experiments (linearly spaced)
trainX = linspace(0,1, nTrainSamples)';
trainY = testFunc(trainX);

testX = linspace(0, 1, 200)';

% GP parameters
ell = 0.05;                                                                 % Characteristic length scale of covariance function
sf = 3;                                                                     
sn = 0.2;                                                                   % "Nugget" term, creates additive noise in the GP.

%% Set Distribution Sequence
sequenceMu = linspace(0.3, 0.7, 20)';
sequenceStd = linspace(0.1, 0.01, 20);
[~, nDist] = size(sequenceMu);
%% Gaussian Process of test

% Set GP Model params
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0; 0];
covfunc = {@covMaterniso, 3}; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.2; hyp.lik = log(sn);

% Optimize the GP params
 hyperParam = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                 covfunc, likfunc, trainX, trainY); 
             
[testY, testYstd]= gp(...
        hyperParam, @infGaussLik, meanfunc, covfunc, likfunc, trainX,...
        trainY, testX...
        );

% Make video
v = VideoWriter('sequentialSurrogate_smooth_2.avi');
open(v);


figure(1)
subplot(3, 1, 1)
plot(testX, testFunc(testX), 'k', 'LineWidth', 2.0, 'LineStyle', '-.')
hold on
scatter(trainX, trainY, 20, 'b', 'filled')

dy = 3*testYstd;
y = testY;
x = testX;

h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 .9],'linestyle','none');
set(h,'facealpha',.5)
line(x,y)
xlabel("X")
ylabel('f(X)')
title('Initial Gaussian Process')
%% Sample Possible Future Gaussian Process Models
nSamp = 100;
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, meanfunc, covfunc, likfunc, ...
    xx, yy, zz);


%% Perform Active Learning
ss = sobolset(1);
ps = scramble(ss,'MatousekAffineOwen');

inputSampler1 = @(numSamp) normrnd(0.6, 0.05, [numSamp, 1]);

gpmlWrap1 = ActiveLearning(gpHandle, trainX,...
                trainY, testX, inputSampler1);
           

XXXX = linspace(0, 1, 100)';
newXXXX = normrnd(0.6, 0.05, [1000, 1]);
numAddSamples = 10;
numTestSamples = 1000;
tic
for iDistribution = 1:nDist
    subplot(3, 1, 2)
    area(XXXX, normpdf(XXXX, sequenceMu(iDistribution), sequenceStd(iDistribution)),'FaceColor','flat')
    hold
    xlim([0,1])
    title("Filter Distribution")
    ylabel('PDF')
    hold off
    inputSampler1 = @(numSamp) normrnd(sequenceMu(iDistribution), sequenceStd(iDistribution), [numSamp, 1]);
    gpmlWrap1.inputSampler = inputSampler1;
    for iIteration = 1:numAddSamples
        subplot(3, 1, 3)
        tic
        [selection1, stdSet1, ~] = gpmlWrap1.activeLearn(XXXX, numTestSamples);
        toc
        yy1 = testFunc(XXXX(selection1));
        gpmlWrap1 = gpmlWrap1.updateTrainingData(XXXX(selection1), yy1);
        [mu1, std1] = gpmlWrap1.predict(XXXX);
        niceErrorbar(XXXX, mu1, std1)
        hold on
        title("Dynamic Gaussian Process")
        plot(XXXX, mu1, 'k', 'LineWidth', 2.0, 'LineStyle', '-.')
        scatter(gpmlWrap1.trainInput, gpmlWrap1.trainOutput, 20, 'b', 'filled')
        ylabel("f(X)")
        xlabel("X")
        hold off
        frame = getframe(gcf);
        for ii = 1:5
            writeVideo(v,frame);
        end
        
        [selection1, stdSet1, ~] = gpmlWrap1.activeUnlearn(numTestSamples);

        gpmlWrap1 = gpmlWrap1.removeTrainingData(selection1);
        [mu1, std1] = gpmlWrap1.predict(XXXX);
        niceErrorbar(XXXX, mu1, std1)
        hold on
        title("Dynamic Gaussian Process")
        plot(XXXX, mu1, 'k', 'LineWidth', 2.0, 'LineStyle', '-.')
        scatter(gpmlWrap1.trainInput, gpmlWrap1.trainOutput, 20, 'b', 'filled')
        ylabel("f(X)")
        xlabel("X")
        xlim([0, 1])
        hold off
        frame = getframe(gcf);
        for ii = 1:5
            writeVideo(v,frame);
        end
    end
end
toc

close(v);




function y = testFunc(x)
    y = sin(x*20);
end

