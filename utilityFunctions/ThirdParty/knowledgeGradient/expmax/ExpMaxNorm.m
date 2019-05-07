function value = ExpMaxNorm(mu,sigma)
% This function calculates \Exp[\max_x \mu_x + \sigma_x Z_x] where Z is a
% standard multivariate normal random variable using numerical integration.
% See my research notes on 3/14/2007 pages 1 and 2 for a derivation of this
% algorithm.

% The algorithm treats two cases.  First case is when all the sigma_x are
% strictly positive.  Second case is when at least one sigma_x is zero.
zeroSigma = find(sigma==0);
if (isempty(zeroSigma))
    f = @(z) 1-MaxNormCdf(z,mu,sigma)-MaxNormCdf(-z,mu,sigma);
    A = 0;
else
    nonzeroSigma = find(sigma>0);
    A = max(mu(zeroSigma));
    f = @(z) 1-MaxNormCdf(z,mu(nonzeroSigma)-A,sigma(nonzeroSigma));
end

value = A + quadgk(f,0,Inf);
