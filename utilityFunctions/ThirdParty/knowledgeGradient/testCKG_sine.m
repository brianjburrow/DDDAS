clear all
close all
clc

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end
start_number = 30;
test_number  = 1000;

xTrain = linspace(-6*pi(), 6*pi(), start_number);
yTrain = -sin(xTrain);

xTest = linspace(-6*pi(), 6*pi(), test_number);
yTest = -sin(xTest);

meanfunc = [];
covfunc = {@covMaterniso, 3}; ell = 1; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 1; hyp.lik = log(sn);


hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                covfunc, likfunc, xTrain', yTrain'); % optimized hyperparams

hyp3 = hyp2;

[mu s2] = gp(...
    hyp3, @infGaussLik, meanfunc, covfunc, likfunc, xTrain',...
    yTrain', xTest'...
    );



 f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
 
xlim([0,1])
subplot(2, 1, 1)
errorbar(xTest, mu, 2*s2)
hold on
plot(xTrain, yTrain, '+')

Sigma = diag(s2);

mu = mu';


N = 10;
noisevar = 0.001 * ones(1, length(s2));
[xkg,maxLogQ,logQ] = CorrelatedNormalKGWithSigma(mu, Sigma, noisevar, N);
hold on
subplot(2, 1, 2)

disp(xTest(xkg))

plot(xTest, logQ./abs(maxLogQ))
hold on
scatter(xTest(xkg), maxLogQ/abs(maxLogQ))
close all
counter = 1;
clear gp

trueMax = xTest(Argmax(yTest));
yMax = max(yTest);
for idx = 1:5
    disp(xTrain)
    ykg = -sin(xTest(xkg));
    xTrain = [xTrain, xTest(xkg)];
    yTrain = [yTrain, ykg];
    [mu s2] = gp(...
    hyp3, @infGaussLik, meanfunc, covfunc, likfunc, xTrain',...
    yTrain', xTest'...
    );
    Sigma = diag(s2);
    mu = mu';
    [xkg,maxLogQ,logQ] = CorrelatedNormalKGWithSigma(mu, Sigma, noisevar, N);
    subplot(5, 2, counter)
    errorbar(xTest, mu, 3*s2)
    hold on
    plot(xTrain, yTrain, '+')
    plot(xTest, yTest)
    scatter(xTest(xkg), ykg, 'o')
    scatter(trueMax, yMax, 'x')
    counter = counter + 1;
    subplot(5, 2, counter)
    plot(xTest, logQ./abs(maxLogQ))
    hold on
    scatter(xTest(xkg), maxLogQ/abs(maxLogQ))
    counter = counter + 1;
    clear gp
end


