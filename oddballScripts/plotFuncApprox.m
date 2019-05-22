clear all
close all
clc

subplot(3, 1, 1)
xx = linspace(0, 10, 1000);
yy = fun(xx);


plot(xx, yy, 'blue','LineWidth', 2.0)
hold on
xlim([0, 10])
ylim([0, 10])
ylabel('\epsilon')
legend('Truth')

subplot(3, 1, 2)
plot(xx, yy, 'blue', 'LineWidth', 2.0)
hold on
libraryIn = linspace(0, 10, 45)';
libraryOut = fun(libraryIn);

% Set GP mean, covariance, and likelihood functions
ell = 0.05;                                                                 % Characteristic length scale of covariance function
sf = 3;                                                                     
sn = 0.2; 
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0; 0];
covfunc = {@covMaterniso, 3}; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.2; hyp.lik = log(sn);

% Perform optimization

hyperParam = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, ...
                 covfunc, likfunc, libraryIn, libraryOut); 

%% Set up instance of the dynamicGP class
% set up GPML function handle
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, ...
    meanfunc, covfunc,...
    likfunc, xx,...
    yy, zz);

[yyy, yerr] = gpHandle(libraryIn, libraryOut, xx');

niceErrorbar(xx', yyy, 4*yerr)
hold on
plot(xx', yyy, 'k', 'LineWidth', 2.0, 'LineStyle', '-.')
xlim([0, 10])
ylim([0, 10])
plot(xx, yy, 'blue', 'LineWidth', 1.0)
ylabel('\epsilon')
f=get(gca,'Children');
legend([f(1),f(2), f(3)],'Truth', 'Approx. 1', 'Est. Error')

%% subplot 3
subplot(3, 1,3)

libraryIn = linspace(0, 10, 10)';
libraryOut = fun(libraryIn);
gpHandle = @(xx, yy, zz) gp(hyperParam, @infGaussLik, ...
    meanfunc, covfunc,...
    likfunc, xx,...
    yy, zz);

[yyy, yerr] = gpHandle(libraryIn, libraryOut, xx');

niceErrorbar(xx', yyy, 4*yerr)
hold on
plot(xx', yyy, 'k', 'LineWidth', 2.0, 'LineStyle', ':')
plot(xx, yy, 'blue', 'LineWidth', 1.0)
xlim([0, 10])
ylim([0, 10])
ylabel('\epsilon')
xlabel('d')
f=get(gca,'Children');
legend([f(1),f(2), f(3)],'Truth', 'Approx.n 2', 'Est. Error')

function y = fun(x)
    y = x + 0.1*x.*sin(5*x);
end