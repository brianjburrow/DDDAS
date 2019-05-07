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
function Sigma = PowerExponentialCovarianceMatrix(cov0, scale, p, MPerDim)
d = length(scale);
M = MPerDim^d;
Sigma = zeros(M); % Preallocate
for i=1:M
    izd = ZToZd(i,MPerDim,d);
    for j=1:M
	jzd = ZToZd(j,MPerDim,d);
	dist = abs(izd-jzd);
	Sigma(i,j) = cov0 * exp(-sum((dist./scale).^p));
    end
end
