function [x,y,variohat,vario]=TestEmpiricalVariogram(cov0,scale,M)
x=[1:M];
p=2;
Sigma=PowerExponentialCovarianceMatrix(cov0, scale, p, M);
mu = zeros(size(x));
y=mvnrnd(mu,Sigma);
variohat = EmpiricalVariogram(x,y);
d = [0:M-1];
vario = TrueVariogram(cov0,scale,M);
[MLEcov0,MLEscale]=EstimateHyperparametersNoNoise(x',y');
MLEvario = TrueVariogram(MLEcov0,MLEscale,M);
plot(d,variohat,d,vario,d,MLEvario);
xlabel('lag')
ylabel('variogram(lag)')
legend('empirical','true','MLE');
title(sprintf('TestEmpricialVariogram: truth(cov0=%g scale=%g) MLE(cov0=%g scale=%g) M=%d', cov0, scale, MLEcov0, MLEscale, M));
