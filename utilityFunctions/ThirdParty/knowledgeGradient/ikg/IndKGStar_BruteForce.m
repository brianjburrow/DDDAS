% Implements the extended independent normal knowledge gradient policy, also called
% concave-enveloped.  Does not use a diffusion approximation, but instead
% explicitly calculates the KG value of measuring an integral number of times,
% for integers between 1 and nleft.
% mu should be an M-vector,
% beliefvar should be an M vector, containing belief variances.
% noisevar should either be a scalar homogeneous noise variance, or an M-vector of the noise variances.
% nleft should be an integer.
function [xkg,maxLogSlope,n]=ExtendedIndKG(mu, beliefvar, noisevar, nleft)
	maxNleft = 10;
	nleft = min(nleft,maxNleft);
	inner_xkg = zeros(1,nleft);
	maxLogQ = zeros(1,nleft);
	for n=1:nleft
		[inner_xkg(n),maxLogQ(n)] = IndKG(mu,beliefvar,noisevar/n);
	end
	logSlope = maxLogQ - log(1:nleft);
	%plot(logSlope)
	[maxLogSlope,bestn] = max(logSlope);
	xkg = inner_xkg(bestn);
end
