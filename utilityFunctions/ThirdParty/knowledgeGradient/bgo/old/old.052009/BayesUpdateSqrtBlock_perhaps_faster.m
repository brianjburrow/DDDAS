% [newmu, newSigmaSqrt] = BayesUpdateSqrtBlock(mu, SigmaSqrt, x, y, noisevar)
% The inputs x, y are vectors of the same length.
% The input noisevar can either be a scalar or the same length as x and y.
% mu is an M-vector.
% SigmaSqrt is an MxM matrix, lower triangular.  Possibly sparse?
% Does not do any post-processing to fix perfect measurements
function [newmu, newSigmaSqrt] = BayesUpdateSqrtBlock(mu, SigmaSqrt, x, y, noisevar)
assert(size(SigmaSqrt,1)==size(SigmaSqrt,2));
assert(size(mu,1)==size(SigmaSqrt,1))
assert(size(mu,2)==1)
n = length(x);
if (length(noisevar)==1)
	noisevar = noisevar*ones(size(x));
end
assert(length(y)==n);
assert(length(noisevar)==n);
newmu = mu;
newSigmaSqrt = SigmaSqrt;
for i=1:n
	%[newmu,newSigmaSqrt] = BayesUpdateSqrt(newmu,newSigmaSqrt,x(i),y(i),noisevar(i));
	tx = x(i);
	ty = y(i);
	tnoisevar = noisevar(i);

	if (tnoisevar == 0 && newSigmaSqrt(tx,tx)==0)
		continue;
	end

	M = length(newmu);
	assert(size(newmu,2)==1);
	F = newSigmaSqrt(tx,:)';
	a = 1/(F'*F + tnoisevar);
	gamma = 1/(1+sqrt(a*tnoisevar));
	K = a*newSigmaSqrt*F;
	residual = ty - newmu(tx);
	newmu = newmu + K*residual;
	newSigmaSqrt = newSigmaSqrt - gamma*K*F';
end
