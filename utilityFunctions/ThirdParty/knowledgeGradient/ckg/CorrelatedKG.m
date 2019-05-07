% This function is currently under construction.  It is intended to implement
% the Correlated KG method, taking either a covariance matrix or its square
% root as input.

% Implements the knowledge gradient policy from 
%   Frazier, P., W.B. Powell, and S. Dayanik, "Knowledge-Gradient Methods for
%   Ranking and Selection with Correlated Normal Beliefs." INFORMS Journal on
%   Computing. 2009.
%
% This is the knowledge gradient policy for the ranking and selection problem
% with correlated normal beliefs, and uncorrelated normal measurement noise
% with known variance.  It can also be used for Bayesian global optimization
% with a gridded domain.  
%
% Allows optionally passing the square root of the covariance matrix instead of
% the covariance matrix itself.  If the measurement noise is small or 0,
% updating the beliefs after measurements using the square-root of the
% covariance matrix can reduce errors due to numerical imprecision.

% options: number of measurements to take, error_tol, d
% PF: need to check for row or column vectorness in mu.

% mu should be an M-vector,
% SigmaSqrt should be a lower diagonal cholesky decomposition of Sigma.
% noisevar can either be an M-vector of the variances, or a scalar.  It doesn't
% matter whether it's a row or column vector.
function [xkg,maxLogQ,logQ]=CorrelatedKG(mu, Sigma, noisevar, use_sqrt, do_tiebreaking)
    % PF: make things optional.  How?
    do_tiebreaking = 0;
    use_sqrt = 1;
    verbose = 0;
    if (use_sqrt)
        SigmaSqrt = Sigma; % Pass in SigmaSqrt in place of Sigma.
        clear Sigma;
    end

    if (~use_sqrt && do_tiebreaking)
        warning('Tie-breaking only supported with SigmaSqrt.');
        do_tiebreaking = 0;
    end
    
    M = length(mu);
    if (verbose)
	disp(sprintf('mu=%s', mat2str(mu,3)));
	% PF: what's a cheap way to get the whole diagonal of Sigma from the sqrt?
	%diag_Sigma = use_sqrt ? diag(
	%disp(sprintf('diag(Sigma)=%s', mat2str(diag(Sigma),3)));
    end
    for x=[1:M]
	if (use_sqrt)
		Sigma_x = SigmaSqrt*SigmaSqrt(x,:)';
	else
		Sigma_x = Sigma(x,:);
	end
	Sigma_xx = Sigma_x(x);
	if (length(noisevar)==1)
		noisevar_x = noisevar;
	else
		noisevar_x = noisevar(x);
	end
	denominator = sqrt(Sigma_xx+noisevar_x);
	if (denominator == 0)
		logQ(x) = -Inf;
		if (verbose > 2)
			disp(sprintf('x=%d logQ=-Inf', x));
		end
	else
		sigmatilde = Sigma_x / denominator;
		logQ(x) = LogEmaxAffine(mu,sigmatilde);
		if (verbose > 2)
			disp(sprintf('x=%d logQ=%f sigmatilde=%s', x, logQ(x), mat2str(sigmatilde,3)));
		end
	end
    end

    % If computers had infinite precision, then we could simply examine the
    % logQ vector to see which alternative has the largest KG-factor.
    % Unfortunately, the real (i.e., infinite precision) differences between
    % Q-factors can be very small, making it hard to tell whether a difference
    % in Q-factor is due to a real difference or numerical noise.  If
    % tie-breaking is enabled then this calls tie-breaking code.  Otherwise, an
    % alternative is selected at random from the tied set.  Note that
    % sometimes, especially on the first measurement with a homogeneous prior,
    % there are real ties not due to numerical noise.

    error_tol = 1e-5;
    xkg = Argmax(logQ);
    maxLogQ=max(logQ);
    tied = find(maxLogQ - logQ < error_tol);
    if (do_tiebreaking && use_sqrt && length(tied) > 1)
	[xkg,maxLogQ,logQ]=tiebreak(mu,SigmaSqrt,noisevar,tied);
    elseif (length(tied)>1)
	disp(sprintf('%d alternatives tied for the largest KG factor at error_tol=%g.  The KG decision returned may be incorrect due to numerical imprecision.  Try enabling tie-breaking code.', length(tied),error_tol));
    end
end

% This function breaks ties between alternatives by calculating differences of
% Q-factors, rather than the Q-factors themselves.  This function was created
% at a time when the numerical precision of LogEmaxAffine was lower and hence
% there were more ties due to numerical noise.  Since then, the numerical
% precision has been improved, and there are fewer ties, so this function may
% be overly cautious, and overly sophisticated for what it needs to do.  This
% function currently only supports passing SigmaSqrt, but it could be extended.
function [xkg,maxLogQ,logQ]=tiebreak(mu, SigmaSqrt, noisevar, tied)
    while (length(tied)>1)
	disp(sprintf('%d alternatives tied for largest KG factor at error_tol=%g: %s', length(ntied), error_tol, mat2str(tied)));
	Sigma_xx = sum(SigmaSqrt(tied(1),:).^2);
	base_sigmatilde = SigmaSqrt*SigmaSqrt(tied(1),:)' / sqrt(Sigma_xx+noisevar(tied(1)));
	clear deltasgn logdelta
	deltasgn(1)=0;
	logdelta(1)=-Inf;
	for k=2:length(tied)
		Sigma_xx = sum(SigmaSqrt(tied(k),:).^2);
		base_sigmatilde = SigmaSqrt*SigmaSqrt(tied(k),:)' / sqrt(Sigma_xx+noisevar(tied(k)));
		[logdelta(k),deltasgn(k)] = LogEmaxAffineDiff(mu,base_sigmatilde,mu,sigmatilde);
	end
	worse = find(deltasgn<0);
	better = find(deltasgn>0);
	same = find(deltasgn==0);
	if (length(worse)>1)
		disp(sprintf('alternatives %s are worse than alternative %d, \nwith log(diff)=%s',...
				mat2str(tied(worse)), tied(1), mat2str(logdelta(worse))));
	end
	if (length(better)>1)
		disp(sprintf('alternatives %s are better than alternative %d, \nwith log(diff)=%s',...
				mat2str(tied(better)), tied(1), mat2str(logdelta(better))));
	end
	if (length(same)>1)
		disp(sprintf('alternatives %s are the same as alternative %d',...
				mat2str(tied(same)), tied(1)));
	end
	if (length(better)>0)
		best = max(logdelta(better));
		stilltied = tied(find(best - logdelta(better) < .1));
		xkg = tied(better(Argmax(logdelta(better))));
	else
		stilltied = tied(same);
		xkg = tied(ChooseRandomElement(same));
	end
	if (stilltied(1)==tied(1))
		if (length(stilltied)>1)
			disp(sprintf('After final tie-breaking, we still have %d tied for best', length(stilltied)));
			disp(sprintf('tied: %s', mat2str(stilltied)));
		end
		break; % and return
	else
		tied=stilltied; % Repeat
	end
    end
end
