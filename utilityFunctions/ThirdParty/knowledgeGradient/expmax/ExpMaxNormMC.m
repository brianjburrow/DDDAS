function [value, err] = ExpMaxNormMC(numpoints, mu,sigma)
% This function calculates \Exp[\max_x \mu_x + \sigma_x Z_x] where Z is a
% standard multivariate normal random variable using Monte Carlo.  Instead
% of using this function, use ExpMaxNorm(mu,sigma), which uses numerical
% integration and is faster for the same accuracy.

% This function only gives a Monte Carlo estimate, and the error is also 
% returned.  We use the qsimvnef function,
% which I downloaded off the internet, and calculates \Exp[f(Z)] for a
% multivariate normal random variable with a covariance matrix that you can
% specify.  numpoints tells how many points to use in the estimate.  The more
% points the more accurate but slower running is the function. 

% M is the dimension of the multivariate normal, a and b are lower and upper
% bounds on the integration, C is the covariance matrix.
M = length(mu);
a = -inf*ones(1,M);
b = inf*ones(1,M);
C = eye(M);
f = @(x) max(mu + sigma.*x');
[ p, e, value, err ] = qsimvnef( numpoints, C, a, b, f );
% p is the probability that the multivariate normal is in the bounds given by a
% and b.  Thus, p will be 1.  e is the error on the estimate of p.  We need
% neither p nor e.  Instead, the third argument is the expectation we are
% interested in, and the fourth argument is its error.
