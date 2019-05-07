clear all
close all
clc
%% Set up test parameters
nTrainSamples = 14;                                                         % Original size of design of experiments (linearly spaced)
trainX = linspace(0,1, nTrainSamples)';
trainY = testFunc(trainX);

testX = linspace(0, 1, 200)';

% GP parameters
ell = 0.05;                                                                 % Characteristic length scale of covariance function
sf = 3;                                                                     
sn = 0.2;                                                                   % "Nugget" term, creates additive noise in the GP.
%% Gaussian Process of test

% Set GP Model params
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0; 0];
covfunc = {@covMaterniso, 3}; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.2; hyp.lik = log(sn);

% Optimize the GP params
 hyperParam = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                 covfunc, likfunc, trainX, trainY); 

% Optimize the GP params
% [logLike, hyperParam, postStruct] = gp(...
%         hyp, @infGaussLik, meanfunc, covfunc, likfunc, trainX,...
%         trainY);
    
[testY, testYstd]= gp(...
        hyperParam, @infGaussLik, meanfunc, covfunc, likfunc, trainX,...
        trainY, testX...
        );


figure(1)
subplot(3,2,[1,2])
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
title('Starting GP')
%% Sample Possible Future Gaussian Process Models
nSamp = 100;
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, meanfunc, covfunc, likfunc, ...
    xx, yy, zz);
%[altGPset, nAlt] = sampleAppendedGP(gpHandle, testX, trainX, trainY, testX, nSamp);

%testSampleAppendedGP(altGPset, testX, nAlt, nSamp)

%% Perform Active Learning
ss = sobolset(1);
ps = scramble(ss,'MatousekAffineOwen');
inputSampler1 = @(numSamp) net(ps, numSamp);
inputSampler2 = @(numSamp) normrnd(0.6, 0.05, [numSamp, 1]);

gpmlWrap1 = ActiveLearning(gpHandle, trainX,...
                trainY, testX, inputSampler1);
            
gpmlWrap2 = ActiveLearning(gpHandle, trainX,...
                trainY, testX, inputSampler2);

XXXX = linspace(0, 1)';
newXXXX = normrnd(0.6, 0.05, [1000, 1]);
numAddSamples = 20;
numTestSamples = 1000;
tic

for iIteration = 1:numAddSamples
    subplot(3,2,3)
    [mu1, std1] = gpmlWrap1.predict(XXXX);
    niceErrorbar(XXXX, mu1, std1)
    hold on
    title("Uniform Prior")
    ylabel("f(X)")
    xlabel('X')
    
    subplot(3,2,4)
    [mu2, std2] = gpmlWrap2.predict(XXXX);
    niceErrorbar(XXXX, mu2, std2)
    hold on
    title("Gaussian Prior")
    ylabel("f(X)")
    xlabel('X')
    
    subplot(3,2,5)

    [selection1, stdSet1, ~] = gpmlWrap1.activeLearn(XXXX, numTestSamples);
    yy1 = testFunc(XXXX(selection1));
    plot(XXXX, stdSet1)
    hold on
    ylabel("Avg Std Reduction")
    xlabel("Alternative X")
    hold off

    subplot(3, 2, 6)
    newXXXX = normrnd(0.6, 0.05, [1000, 1]);
    [selection2, stdSet2, ~] = gpmlWrap2.activeLearn(XXXX, numTestSamples);
    plot(XXXX, stdSet2);
    hold on
    ylabel("Avg Std Reduction")
    xlabel("Alternative X")
    hold off
    yy2 = testFunc(XXXX(selection2));

    
    subplot(3,2,3)
    scatter(XXXX(selection1), yy1, 'filled', 'g')
    hold off
    
    subplot(3,2,4)
    scatter(XXXX(selection2), yy2, 'filled', 'g')
    hold off
    %% update GPs
    gpmlWrap1.updateTrainingData(XXXX(selection1), yy1);
    gpmlWrap2.updateTrainingData(XXXX(selection2), yy2);
    pause(0.05)
end
toc





XXXX = linspace(0, 1)';
newXXXX = normrnd(0.6, 0.05, [1000, 1]);
numTestSamples = 1000;
tic
gpmlWrap1.sortTrainingData();
gpmlWrap2.sortTrainingData();
numRemoveSamples = 10;
for iIteration = 1:numRemoveSamples
    subplot(3,2,3)
    [mu1, std1] = gpmlWrap1.predict(XXXX);
    niceErrorbar(XXXX, mu1, std1)
    hold on
    title("Uniform Prior")
    ylabel("f(X)")
    xlabel('X')
    
    subplot(3,2,4)
    [mu2, std2] = gpmlWrap2.predict(XXXX);
    niceErrorbar(XXXX, mu2, std2)
    hold on
    title("Gaussian Prior")
    ylabel("f(X)")
    xlabel('X')
    
    subplot(3,2,5)

    [selection1, stdSet1, ~] = gpmlWrap1.activeUnlearn(numTestSamples);
    yy1 = testFunc(XXXX(selection1));
    plot(gpmlWrap1.trainInput, stdSet1)
    hold on
    ylabel("Avg Std Reduction")
    xlabel("Alternative X")
    xlim([0, 1])
    hold off

    subplot(3, 2, 6)
    newXXXX = normrnd(0.6, 0.05, [1000, 1]);
    [selection2, stdSet2, ~] = gpmlWrap2.activeUnlearn(numTestSamples);
    plot(gpmlWrap2.trainInput, stdSet2);
    hold on
    ylabel("Avg Std Reduction")
    xlabel("Alternative X")
    xlim([0,1])
    hold off
    yy2 = testFunc(XXXX(selection2));

    
    subplot(3,2,3)
    scatter(gpmlWrap1.trainInput(selection1), yy1, 'filled')
    hold off
    
    subplot(3,2,4)
    scatter(gpmlWrap2.trainInput(selection2), yy2, 'filled')
    hold off
    %% update GPs
    gpmlWrap1.removeTrainingData(selection1);
    gpmlWrap2.removeTrainingData(selection2);
    pause(0.5)
end
toc




aoeue




function y = testFunc(x)
    y = sin(x*20);
end

