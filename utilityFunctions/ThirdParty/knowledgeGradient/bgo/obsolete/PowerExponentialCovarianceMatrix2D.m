%Sigma = PowerExponentialCovarianceMatrix2d(cov0, scale, p, MPerDim)
%
%cov0 is a scalar
%scale is a 2-dimensional vector (also called 1./theta in the dace literature).
%p is a scalar
%MPerDim is a scalar
function Sigma = PowerExponentialCovarianceMatrix2d(cov0, scale, p, MPerDim)
% warning('This function is deprecated in favor of PowExpCov, which is more accurate.');
assert(length(scale)==2);
M = MPerDim^2;
Sigma = zeros(M); % Preallocate
for i=1:M
    [i1zd i2zd] = ind2sub([MPerDim MPerDim],i);
    [j1zd j2zd] = ind2sub([MPerDim MPerDim],1:M);
    dist1 = abs(i1zd-j1zd)/scale(1);
    dist2 = abs(i2zd-j2zd)/scale(2);
    Sigma(i,:) = cov0 * exp(-(dist1.^p + dist2.^p));
end
