% DO NOT USE THIS FUNCTION.  It is better to either use BayesUpdateBlock to
% update a block of measurements as a batch, or to work with the square root of
% the covariance matrix either recursively (BayesUpdateSqrt) or in batch
% (BayesUpdateSqrtBatch). 
%
% [newmu, newSigma] = BayesUpdate(mu, Sigma, x, y, noises)
% The inputs x, y, and noises are scalars.
% Corrects newmu(x) to y and newSigma(x) to 0 on perfect measurements.
function [newmu, newSigma] = BayesUpdate(mu, Sigma, x, y, noises)
	[newmu,newSigma] = BayesUpdateBlock(mu,Sigma,x,y,noises);
end

