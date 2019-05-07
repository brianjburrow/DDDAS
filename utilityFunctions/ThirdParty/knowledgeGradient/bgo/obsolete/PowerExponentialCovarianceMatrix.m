%Sigma = PowerExponentialCovarianceMatrix(cov0, scale, p, MPerDim)
%
%cov0 is a scalar
%scale is a d-dimensional vector (also called 1./theta in the dace literature).
%p is a scalar
%MPerDim is a scalar
%
%I tested this code just by running,
%PowerExponentialCovarianceMatrix(.5,[.5,.5],2,3)
%and making sure it didn't report any errors, and eyeballing the values
%returned.
% This function is deprecated in favor of PowExpCov, which is more accurate.
function Sigma = PowerExponentialCovarianceMatrix(cov0, scale, p, MPerDim)
% We pass off to different routines for speed.
d = length(scale);
if (d==1)
	Sigma = PowerExponentialCovarianceMatrix1D(cov0,scale,p,MPerDim);
elseif (d==2)
	Sigma = PowerExponentialCovarianceMatrix2D(cov0,scale,p,MPerDim);
elseif (d==3)
	Sigma = PowerExponentialCovarianceMatrix3D(cov0,scale,p,MPerDim);
else
	Sigma = PowerExponentialCovarianceMatrixAnyD(cov0,scale,p,MPerDim);
end
