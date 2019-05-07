%{
Creates a finite dimensional covariance matrix for a gaussian random function
with a power exponential covariance function on a line.  This covariance
function is Cov(x,y) = cov0 * R(x-y), where R is the correlation function
defined by R(d) = exp(-|d/scale|^p).  This computation gives the infinite
dimensional distribution, from which the function computes the finite
dimensional distribution for the points {0,\delta,...,1-\delta,1}, where
\delta=1/(M-1).  cov0 is called beta and scale is called alpha in the paper.
%}
function Sigma = PowerExponentialCovarianceMatrix1D(cov0, scale, p, M, nowarn)
% warning('This function is deprecated in favor of PowExpCov, which is more accurate.');
assert(length(scale)==1);
Sigma = zeros(M); % Preallocate
for i=[1:M]
    j = 1:M;
    d = abs(i-j); % Before this was delta*abs(i-j), where delta=1/M-1.  Instead, multiply scale by M-1.
    Sigma(i,j) = cov0 * exp(-(d/scale).^p);
end
