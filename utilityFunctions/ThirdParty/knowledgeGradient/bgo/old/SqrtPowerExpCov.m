function SigmaSqrt = SqrtPowerExpCov(cov0, scale, p, MPerDim)
% We pass off to different routines for speed.
d = length(scale);
if (d==1)
	SigmaSqrt = SqrtPowerExpCov1D(cov0,scale,p,MPerDim);
elseif (d==2)
	SigmaSqrt = SqrtPowerExpCov2D(cov0,scale,p,MPerDim);
else
	error('SqrtPowerExpCov only currently supports dimensions 1 and 2');
end

% Test the quality of the cholesky decomposition.
err = max(max(abs(PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim)-SigmaSqrt*SigmaSqrt')));
if (err > 1e-2)
	warning(sprintf('SqrtPowerExpCov: err=max(Sigma-SigmaSqrt^2)=%g, err/cov0=%g',err,err/cov0));
end
