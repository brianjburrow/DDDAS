% Given a prior multivariate normal prior with mean vector mu and covariance
% matrix Sigma, as well as a measurement y of entry x with measurement noise
% variance noisevar, provides the posterior.  The inputs x, y, and noisevar are
% scalars.  Corrects newmu(x) to y and newSigma(x) to 0 on perfect
% measurements (measurements with noisevar=0).
%
% Use this function at your own risk.  Often it is better to either use
% BayesUpdateBlock to update a block of measurements as a batch, or to work
% with the square root of the covariance matrix either recursively
% (BayesUpdateSqrt) or in batch (BayesUpdateSqrtBatch). 
%
% [newmu, newSigma] = BayesUpdate(mu, Sigma, x, y, noisevar)
function [newmu, newSigma] = BayesUpdate(mu, Sigma, x, y, noisevar)
	[newmu,newSigma] = BayesUpdateBlock(mu,Sigma,x,y,noisevar);
end

