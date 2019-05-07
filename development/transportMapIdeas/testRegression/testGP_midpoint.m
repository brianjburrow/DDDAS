clear all
clc
close all

polyOrders = 1;
dim = 1;  % Dimensionality of parameter space
lb = 0;
ub = 500;
numSlice = 500;
numSamples = 200;
noiseVar = 10;
%samps = genPosteriorSamples(4);


%% Run Script
XXX = linspace(lb, ub, numSlice*4);
xxx = linspace(lb, ub, numSlice);
testX = linspace(lb, ub, numSlice*8);
samples = zeros(numSamples, numSlice);
figure(1)
subplot(4,1,1)


funct = @(xx) xx*sin(xx/(100 + xx));
for idx = 1:length(xxx)
    yyy(idx) = funct(xxx(idx)) + normrnd(0, noiseVar^0.5);
end

fullGP = fitrgp(xxx', yyy'  ,'KernelFunction','squaredexponential',...
   'PredictMethod', 'exact', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
   struct('AcquisitionFunctionName','expected-improvement-plus'));


 sdGP = fitrgp(xxx', yyy'  ,'KernelFunction','squaredexponential',...
    'PredictMethod', 'sd', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
 
 ficGP = fitrgp(xxx', yyy'  ,'KernelFunction','squaredexponential',...
    'PredictMethod', 'fic', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));

[YYtrain, errF] = predict(fullGP, XXX');

muR = fitrsvm(XXX', YYtrain,'KernelFunction','gaussian','KernelScale','auto',...
    'Standardize',true, 'OptimizeHyperparameters','expected-improvement-plus');

SigmaR = fitrsvm(XXX', errF,'KernelFunction','gaussian','KernelScale','auto',...
    'Standardize',true, 'OptimizeHyperparameters','expected-improvement-plus');


figure(1)
subplot(5, 1, 1)
scatter(xxx, yyy)
ylim([-1000, 1000])

subplot(5, 1, 2)
tic
[YY, errF] = predict(fullGP, testX');
F = toc;
errorbar(testX, YY, 3*errF)
ylim([-1000, 1000])

subplot(5, 1, 3)
tic
[YYsd, errSd] = predict(fullGP, testX');
SD = toc;
errorbar(testX, YYsd, 3*errSd)
ylim([-1000, 1000])

subplot(5, 1, 4)
tic
[YYfic, errfic] = predict(fullGP, testX');
FIC = toc;
errorbar(testX, YYfic, 3*errfic)
ylim([-1000, 1000])

tic
est = predict(muR, testX');
err = predict(SigmaR, testX');
SVM = toc;
subplot(5,1,5)
errorbar(testX, est, 3*err)
ylim([-1000, 1000])

disp(F - F)
disp(FIC - F)
disp(SD - F)
disp(SVM - F)

