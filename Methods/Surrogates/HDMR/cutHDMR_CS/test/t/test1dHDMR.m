clear all 
close all
nSamples = 100;
lowBound = -1.5;
upBound = 1.5;
xxx = linspace(lowBound, upBound, 10000);
yyy = testFunc(xxx);

%% Prepare figure
fig = figure(1);
ax1 = subplot(3, 1, 1);
ylabel('f(x)')
hold on
ylim([-20, 30])
xlim([-1.5, 1.5])
ax2 = subplot(3, 1, 2);
ylabel('f(x)')
hold on
ylim([-20, 30])
xlim([-1.5, 1.5])
ax3 = subplot(3, 1, 3);
ylabel('f(x)')
hold on
xlabel('x')
ylim([-20, 30])
xlim([-1.5, 1.5])



%% Test a global sampling
% Generate initial results
x = linspace(lowBound, upBound, nSamples);
cutPoint = 0;
f0 = testFunc(cutPoint) * ones(size(x));
F0 = testFunc(cutPoint) * ones(size(xxx));
f1 = testFunc(x) - f0;
tic
lpc0 = LegendrePolynomialCurve(5, x', f1');
y = lpc0.predict(xxx') + F0';
toc

% plot those results
set(fig, 'currentaxes', ax1)
plot(xxx,yyy, 'k', 'Linewidth', 2.0);
scatter(x, f0, 10, 'filled', 'k');
scatter(x, f1, 10, 'filled', 'b');
approx = plot(xxx', y, 'g');



%% Test local deterministic sampling
lowBound = 0.25;
upBound = 1.0;
cutPoint = 0.4;

% generate second set of results
x = linspace(lowBound, upBound, nSamples);                                         % Test new design
f0 = testFunc(cutPoint) * ones(size(x));
F0 = testFunc(cutPoint) * ones(size(xxx));
f1 = testFunc(x) - f0;
tic
lpc1 = LegendrePolynomialCurve(5, x', f1');
y = lpc1.predict(xxx') + F0';
toc

% plot those results
set(fig, 'currentaxes', ax2)
plot(xxx,yyy, 'k', 'Linewidth', 2.0);
scatter(x, f0, 10, 'filled', 'k');
scatter(x, f1, 10, 'filled', 'b');
plot(xxx', y, 'g');



%% Test local random sampling
%x = normrnd(0.6, 0.2, 1, nSamples);
x = normrndSobol(0.6, 0.43, nSamples, 1)';

f0 = testFunc(cutPoint) * ones(size(x));
F0 = testFunc(cutPoint) * ones(size(xxx));
f1 = testFunc(x) - f0;
tic
lpc2 = LegendrePolynomialCurve(5, x', f1');
toc
tic
y = lpc2.predict(xxx') + F0';
toc

set(fig, 'currentaxes', ax3)
plot(xxx,yyy, 'k', 'Linewidth', 2.0);
scatter(x, f0, 10, 'filled', 'k');
scatter(x, f1, 10, 'filled', 'b');
plot(xxx', y, 'g');

function output = testFunc(xx)
    % Fifth order Hermite Polynomial
    output = xx.^5 - 10*xx.^3 + 15*xx + 5*xx.*sin(20*xx) + 10*ones(size(xx));
end

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end

function points = normrndSobol(MU, COV, nPoints, nDimensions)
    ss = sobolset(nDimensions);
    stdNormPoints = uniformTOgaussian(ss(2:nPoints+1, :));
    points = MU + COV*stdNormPoints*COV';
end