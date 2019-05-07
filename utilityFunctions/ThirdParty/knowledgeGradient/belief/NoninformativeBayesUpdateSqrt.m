% Does a NoninformativeBayesUpdate with the A matrix and the square root of the
% B matrix.  This should be somewhat more numerically stable than vanilla
% NoninformativeBayesUpdate.
%
% Assumes that B0 = SqrtB0*SqrtB0', so you can use SqrtB0=chol(B0,'lower') to
% compute SqrtB0.  Be careful not to use chol(B0), since this will produce a
% SqrtB0 that satisfies B0 = SqrtB0'*SqrtB0, with the transpose in the wrong
% position.
function [mu1,A1,SqrtB1] = NoninformativeBayesUpdateSqrt(mu0,A0,SqrtB0,x0,lambdax,yhat1)
    % Check that sqrtB0 is lower triangular.
    if(~all(all(SqrtB0 == tril(SqrtB0))))
	error('NoninformativeBayesUpdateSqrt should be called with a lower triangular SqrtB0');
    end
    if (SqrtB0(x0,x0) ~= 0)
	B0 = SqrtB0*SqrtB0';
        mu1 = mu0 + (yhat1 - mu0(x0))*B0(:,x0)/B0(x0,x0);
        A1 = A0 + ((A0(x0,x0)+lambdax)*B0 - A0(:,x0)*B0(x0,:) - B0(:,x0)*A0(x0,:))/B0(x0,x0);

	% We want to update B1 according to 
	% B1 = B0 - B0(:,x0)*B0(x0,:)/B0(x0,x0), 
	% but we want to do it with square roots.  This is nearly the same
	% update that is done for Kalman filtering with an informative prior:
	% newSigma = Sigma - Sigma(:,x0)*Sigma(x0,:)/(Sigma(x0,x0)+noisevar).
	% If we set noisevar to 0 and Sigma to B0, then B1 is obtained as
	% newSigma.  Thus, to do this update with the square roots of B0 and
	% B1, we can use BayesUpdateSqrt passing a 0 noise variance.
	% That code takes a mean vector and covariance matrix and updates it,
	% but we just use the covariance matrix part of it, and throw away the
	% update to the mean.
	tmp = zeros(size(mu0)); % this is our fake mean.
	[tmp, SqrtB1] = BayesUpdateSqrt(tmp, SqrtB0, x0, 0, 0);

	% Even though we are using square roots, and this is more numerically
	% stable than just using B0 and B1 directly, the previous step does not always 
	% get the columns and rows corresponding to the covariance with x0 to
	% be perfectly 0.  We know theoretically that it must be, because now
	% the posterior on the value at x0 is informative, so we explicitly set
	% these to 0 here.  
	SqrtB1(x0,:) = 0;
	SqrtB1(:,x0) = 0;

	% Note to the user: If, in calling this routine, you know from theory
	% that other alternatives will also be informative after the update in
	% your problem, explicitly set these to 0 too outside this function in
	% your code. 
    else
	warning('An alternative with an informative marginal prior is being measured');
        mu1 = mu0 + (yhat1 - mu0(x0))*A0(:,x0)/(A0(x0,x0)+lambdax);
        A1 = A0 - A0(:,x0)*A0(x0,:)/(A0(x0,x0) + lambdax);
        SqrtB1 = SqrtB0;
    end
