%[impl,x,y,extras] = SimFixedHyperSKO(N,sample,dim,n0,p,MPerDim,mu0,cov0,scale,noisevar)
%
% Simulates SKO algorithm through N samples, and collects at each sample time
% 1<=n<=N the measurement and implementation decision.  Uses fixed
% hyperparameters.
%
% N -- number of measurements to take.
% sample(x) -- function that, given a (flattened) measurement location x, returns a
%   noisy sample of the function value at that point.
% dim -- number of dimensions in the search domain.
% n0 -- number of measurements to put in the first stage.
% MPerDim -- number of discretizations in each dimension.  A scalar.
% mu0, cov0, scale, noisevar -- hyperparameters to use.  All but scale are scalars.
%
% impl -- vector of implementation decisions
% x -- vector of measurement decisions
% y -- vector of observations
% extras -- Not used.  Returns NaN.
%
function [impl,x,y,extras] = SimFixedHyperSKO_v3(N,sample,dim,n0,p,MPerDim,mu0,cov0,scale,noisevar)
    if (p~=2)
	error('This function is hard-coded to p=2');
    end
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras = NaN;
    toZd = @(z) ZToZd(z,MPerDim,dim);

    mun = mu0*ones(M,1); % mu0 is a scalar
    [tmp,SigmaSqrt]=PowExpCov(cov0,scale,p,MPerDim);
    clear tmp

    if (n0 > 0)
	% Run a first stage using latin hypercube sampling.
	xd(1:n0,:) = DiscreteLHS(MPerDim,dim,n0);
	for n=1:n0
	    x(n) = ZdToZ(xd(n,:),MPerDim);
	    y(n) = sample(x(n));
	    % In this first stage, the implementation decision is the best of the
	    % observations we have made thus far.  In the second stage, we allow
	    % the implementation decision to range over the whole space, being the
	    % best of the current mu.  This is a discontinuity, and should be kept
	    % in mind.
	    if (n<n0)
		impl(n) = x(Argmax(y(1:n)));
		disp(sprintf('n=%d: SKO x=%d (%s)', n, x(n), mat2str(xd(n,:))));
	    end
        end

        [mun,SigmaSqrt]=BayesUpdateSqrtBlock(mun,SigmaSqrt,x(1:n),y(1:n),noisevar);

        impl(n) = Argmax(mun);
        disp(sprintf('n=n0=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g', ...
		n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
		mun(impl(n)), mat2str(scale,3), cov0, noisevar)); 
    end

    % Run the rest of the measurements
    for n=[n0+1:N]
	% mun, SigmaSqrt have been set up by the previous iteration of this loop.
	sposterior = sqrt(full(diag(SigmaSqrt*SigmaSqrt')));
	if (noisevar == 0)
		% If measurements are perfect, we correct for them and use EGO.
		sposterior(x(1:n-1))=0;
    		x(n) = EGO(-mun,sposterior);
	else
		c = 1; % optional argument to SKO
    		x(n) = SKO(-mun,sposterior,sqrt(noisevar),c);
	end
	xd(n,:) = toZd(x(n));
	y(n) = sample(x(n));

	% Do the update & get the implementation decision.
	[mun,SigmaSqrt]=BayesUpdateSqrt(mun,SigmaSqrt,x(n),y(n),noisevar);
	if (noisevar == 0)
		% Correct for perfect measurements
		mun(x(1:n))=y(1:n);
	end
	impl(n) = Argmax(mun);
	disp(sprintf('n=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g', ...
			n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
			mun(impl(n)), mat2str(scale,3), cov0, noisevar)); 
    end
end
