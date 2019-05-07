% [newmu, newSigmaSqrt] = BayesUpdateSqrt(mu, SigmaSqrt, x, y, noisevar)
% The inputs x, y, and noisevar are scalars.
% mu is a column M-vector.
% SigmaSqrt is an MxM matrix, lower triangular.  Possibly sparse?
% Does not do any post-processing to fix perfect measurements
function [newmu, newSigmaSqrt] = BayesUpdateSqrt(mu, SigmaSqrt, x, y, noisevar)
% Preprocess to get rid of measurements for which we already know the result
% perfectly.
if (noisevar == 0 && SigmaSqrt(x,x)==0)
	if (mu(x)~=y)
		warning(sprintf('Perfect measurment of x=%d had y=%g but mu(x)=%g with y-mu(x)=%g', x,y,mu(x),y-mu(x)));
	end
	newmu = mu;
	newSigmaSqrt = SigmaSqrt;
	return
end

M = length(mu);
assert(size(mu,2)==1); % Must be a column vector
F = SigmaSqrt(x,:)'; % In KaBrSc71, F=SigmaSqrt'*C', with C=e(x,M)', which is equivalent.
a = 1/(F'*F + noisevar); % In KaBrSc71, noisevar is called \Theta.
gamma = 1/(1+sqrt(a*noisevar));
K = a*SigmaSqrt*F;
residual = y - mu(x);
newmu = mu + K*residual;
newSigmaSqrt = SigmaSqrt - gamma*K*F';

% KaBrSc71 says that numerical imprecision is a problem as the condition number
% approaches 10^p, where p is the number of base-10 digits of precision in your
% floating point.  For us, calling eps(10) gives 1.78e-15, which tells me that
% we have about 15 base-10 digits of precision.  To be safer, I warn when we
% get to 1e12.
% PF, comment added later: I commented out this check because it doesn't seem
% to have much meaning.  For example, if we do a perfect measurement SigmaSqrt
% will have a 0 on the diagonal and condition number infinity, but we can still
% do computations.  Is there a better way to check the numerical precision of
% this computation?
%{
if (condest(SigmaSqrt)>1e12)
	disp(sprintf('Condition number of SigmaSqrt is %g', condest(SigmaSqrt)));
end
%}
