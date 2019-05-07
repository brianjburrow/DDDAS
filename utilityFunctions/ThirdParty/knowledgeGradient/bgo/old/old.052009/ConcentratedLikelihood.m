% x and y should be column vectors.
function [likelihood, muhat, sigma2hat, psi] = ConcentratedLikelihood(x, y, theta, p)
CheckColVector(x);
CheckColVector(y);
n = length(x); % number of data points
R = zeros(n);
for i=1:n
    for j=1:n
	% dij is the distance between x(i) and x(j).  This code works
	% in both the case that theta and p are vectors with the same length as
	% the dimension of x, and in the case that they are scalars, and in the
	% case that one is a vector and the other a scalar.
	dij = theta'*((x(i)-x(j)).^p);
	R(i,j) = exp(-dij);
    end
end

Rinv = R^-1;
one = ones(n,1); % column vector
muhat = (one'*Rinv*y)/(one'*Rinv*one);
res = y-one*muhat;
sigma2hat = (res'*Rinv*res)/n;
const = (2*pi*exp(1))^(-n/2);
likelihood = const / (sigma2hat^(n/2)*sqrt(det(R)));
psi = det(R)^(1/n) * sigma2hat;
