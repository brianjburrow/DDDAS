%Sigma = PowerExponentialCovarianceMatrix2d(cov0, scale, p, MPerDim)
%
%cov0 is a scalar
%scale is a 3-dimensional vector (also called 1./theta in the dace literature).
%p is a scalar
%MPerDim is a scalar
function Sigma = PowerExponentialCovarianceMatrix2d(cov0, scale, p, MPerDim)
% warning('This function is deprecated in favor of PowExpCov, which is more accurate.');
assert(length(scale)==3);
M = MPerDim^3;
Sigma = zeros(M); % Preallocate
for i=1:M
    [i1zd i2zd i3zd] = ind2sub([MPerDim MPerDim MPerDim],i);
    [j1zd j2zd j3zd] = ind2sub([MPerDim MPerDim MPerDim],1:M);
    dist1 = abs(i1zd-j1zd)/scale(1);
    dist2 = abs(i2zd-j2zd)/scale(2);
    dist3 = abs(i3zd-j3zd)/scale(3);
    Sigma(i,:) = cov0 * exp(-(dist1.^p + dist2.^p + dist3.^p));
end
