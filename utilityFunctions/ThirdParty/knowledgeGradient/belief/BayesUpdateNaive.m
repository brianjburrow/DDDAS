% [newmu, newSigma] = BayesUpdate(mu, Sigma, x, yhat, lambda_x)
% The input lambda_x is a scalar.
function [newmu, newSigma] = BayesUpdateNaive(mu, Sigma, x, yhat, lambda_x)
    warning('this function does not work very well.  It sometimes introduces numerical imprecision.  Use BayesUpdateBlock instead.');
    s = Sigma(:,x); % Column vector
    d = s(x) + lambda_x;
    if (d==0)
        % This is the case where we already knew perfectly the value at x,
        % and we measured it with 0 precision.
	if (nargout > 1)
		newSigma = Sigma;
	end
        if (abs(mu(x)-yhat)>10*eps(yhat));
		disp(sprintf('mu(x)=%f and yhat=%f differ (diff=%g, eps=%g) on a perfect measurement with perfect a priori knowledge', mu(x), yhat, mu(x)-yhat, eps(yhat)));
	end
        newmu = mu;
        return;
    end
    newmu = mu + s*(yhat-mu(x))/d;
    if (nargout > 1)
    	newSigma = Sigma - s*s'/d;

	% This is a fix to handle numerical imprecision that occurs in the case
	% when lambda_x = 0.  Without this fix, we can sometimes get negative
	% values on the diagonal.  A better way might be to keep the square
	% root of the covariance matrix, as discussed in the literature.
	M = length(s);
	if (lambda_x == 0)
		newSigma(x,:) = zeros(1,M);
		newSigma(:,x) = zeros(M,1);
	end
	mineig = min(eig(Sigma));
	loop = 0;
	while (mineig<0 && loop<10)
		newSigma = newSigma + (-mineig+eps(-mineig))*eye(M);
		disp(sprintf('BayesUpdate: fixing positive semi-definiteness of covariance matrix, adding %g, loop %d', mineig, loop));
		mineig = min(eig(newSigma));
		loop = loop + 1;
	end
	if(mineig<0)
		disp(sprintf('Even after fixing positive-semidefiniteness, Sigma has min(eig)=%g', mineig)); 
	end
	assert(all(all(isreal(Sigma))));
    end
