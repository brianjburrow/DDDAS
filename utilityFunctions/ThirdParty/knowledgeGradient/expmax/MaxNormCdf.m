% Computes the probability of the event 
% {max_x mu_x + sigma_x Z_x <= alpha}, where the Z_x are independent
% standard normals.  Assumes that all of the sigma_x are strictly positive.
% The passed mu and sigma should be row vectors.  The passed alpha may be a
% scalar or a row vector.  If the passed alpha is a row vector, it computes
% the tail probability for each element of the vector, and returns these
% tail probabilities in a row vector.
%
% The algorithm used is
% P{max_x mu_x + sigma_x Z_x <= alpha} 
% = P{mu_x + sigma_x Z_x <= alpha for all x}
% = \Product_x P{mu_x + sigma_x Z_x <= alpha} by independence
% = \Product_x P{Z_x <= (alpha - mu_x)/sigma_x} since sigma_x > 0.
function P = MaxNormCdf(alpha,mu,sigma)
for i=[1:length(alpha)]
    z(:,i) = (alpha(i) - mu)./sigma;
end
% The second "1" argument to prod is necessary because, if mu and sigma
% each have only one element, then z will be a row vector, which is the
% degenerate case of an Mxlength(alpha) matrix.  In general, we want prod 
% to multiply along each column of z, but when length(alpha)=1 its default
% behavior is to multiply the entries of the row vector, supplying just one
% number rather than M.  Supplying the second "1" argument enforces
% multiplication along each column of z, even when z is a row vector.
P = prod(normcdf(z),1);

% One test case is 
% >> MaxNormCdf([1,2,3,4,5],[0],[1])-normcdf([1,2,3,4,5])
% This should come out all zeros.
%
% Another test case is looking at MaxNormCdf([1,2,3,4,5,6],mu,sigma)
% for any valid mu, sigma.  This should return a vector increasing to 1.
% >> mu=rand(1,10);
% >> sigma=rand(1,10);
% >> MaxNormCdf([1,2,3,4,5,6],mu,sigma)