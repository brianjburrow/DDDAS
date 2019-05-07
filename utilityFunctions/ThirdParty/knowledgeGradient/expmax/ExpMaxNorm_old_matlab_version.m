function value = ExpMaxNorm(mu,sigma)
% This function calculates \Exp[\max_x \mu_x + \sigma_x Z_x] where Z is a
% standard multivariate normal random variable using numerical integration.
% See my research notes on 3/14/2007 pages 1 and 2 for a derivation of this
% algorithm.

% The implementation of this code only uses the quad function, which does not
% integrate to infinity, rather than the fancier function quadgk that can
% integrate to infinity.

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

% We want to integrate to infinity, but the quad function can't integrate 
% all the way to infinity.  One solution is just to integrate to a big
% number.;
% Integration limit.  This is a hack
limit=10^8;
value = A + quad(f,0,limit);
