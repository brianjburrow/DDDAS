%[impl,x,y,extras] = SimAdaptiveSKO(N,sample,dim,n0a,n0b,MPerDim,MLE)
%
% Simulates SKO algorithm through N samples, and collects at each sample time
% 1<=n<=N the measurement and implementation decision.  Uses the passed 
% MLE to estimate hyperparameters, including the noise variance.
%
% N -- number of measurements to take.
% sample(x) -- function that, given a (flattened) measurement location x, returns a
%   noisy sample of the function value at that point.
% dim -- number of dimensions in the search domain.
% n0a -- number of measurements to put in the first stage latin hypercube.
% n0b -- number of measurements to put in the second stage, when we measure
%	again the n0b best points.
% MPerDim -- number of discretizations in each dimension.  A scalar.
% MLE -- function handle of the form,
%	[cov0hat, noisevarhat, scalehat, betahat] = MLE(xd, y)
%
% impl -- vector of implementation decisions
% x -- vector of measurement decisions
% y -- vector of observations
% extras.cov0 -- array of estimated cov0
% extras.mu0 -- array of estimated mu0
% extras.scale -- array of estimated scale
% extras.noisevar -- array of estimated noise variances
%
function [impl,x,y,extras] = SimAdaptiveSKO(N,sample,dim,n0a,n0b,MPerDim,MLE)
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras.cov0 = NaN*ones(1,N);
    extras.mu0 = NaN*ones(1,N);
    extras.scale = NaN*ones(dim,N);
    extras.noisevar = NaN*ones(1,N);

    toZd = @(z) ZToZd(z,MPerDim,dim);

    % Sample n0a points using latin hypercube sampling.
    n0 = n0a+n0b;
    assert(n0a>0);
    xd(1:n0a,:) = DiscreteLHS(MPerDim,dim,n0a);
    for n=1:n0a
	x(n) = ZdToZ(xd(n,:),MPerDim);
	y(n) = sample(x(n));
	if (n<n0)
		impl(n) = x(Argmax(y(1:n)));
		disp(sprintf('n=%d: SKO x=%d (%s)', n, x(n), mat2str(xd(n,:))));
	end
    end
    % Take the next n0b points from the best points measured thus far.
    if (n0b>0)
	assert(n0b<n0a);
	[besty,besti]=sort(y(1:n0a),2,'descend');
	x(n0a+1:n0)=x(besti(1:n0b));
	for n=n0a+1:n0
        y(n) = sample(x(n));
		xd(n,:) = toZd(x(n));
		if (n<n0)
			impl(n) = x(Argmax(y(1:n)));
			disp(sprintf('n=%d: SKO x=%d (%s)', n, x(n), mat2str(xd(n,:))));
		end
	end
    end

    [extras.cov0(n), extras.noisevar(n), extras.scale(:,n), extras.mu0(n)] = MLE(xd(1:n,:), y(1:n));
    mun = extras.mu0(n)*ones(M,1);,
    assert(length(extras.scale(:,n))==dim);
    [tmp,SigmaSqrt]=PowExpCov(extras.cov0(n),extras.scale(:,n),2,MPerDim);
    [mun,SigmaSqrt]=BayesUpdateSqrtBlock(mun,SigmaSqrt,x(1:n),y(1:n),extras.noisevar(n));

    impl(n) = Argmax(mun);
    disp(sprintf('n=n0=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g mu0=%g', ...
		n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
		mun(impl(n)), mat2str(extras.scale(:,n),3), extras.cov0(n), ...
		extras.noisevar(n), extras.mu0(n))); 

    % Run the rest of the measurements
    for n=[n0+1:N]
	% mun, SigmaSqrt have been set up by the previous iteration of this loop.
	% If measurements are perfect, we should correct for them.  This is not
	% coded yet.
	if(extras.noisevar(n-1) == 0)
		warning('We need to correct measurements when they are perfect');
	end
	sposterior = full(sqrt(diag(SigmaSqrt*SigmaSqrt')));
	c = 1; % optional argument to SKO
	x(n) = SKO(-mun,sposterior,sqrt(extras.noisevar(n-1)),c);

	xd(n,:) = toZd(x(n));
	y(n) = sample(x(n));

	% Do the update & get the implementation decision.
	[extras.cov0(n), extras.noisevar(n), extras.scale(:,n), extras.mu0(n)] = MLE(xd(1:n,:), y(1:n));
	mun = extras.mu0(n)*ones(M,1);
	[tmp,SigmaSqrt]=PowExpCov(extras.cov0(n),extras.scale(:,n),2,MPerDim);
	[mun,SigmaSqrt]=BayesUpdateSqrtBlock(mun,SigmaSqrt,x(1:n),y(1:n),extras.noisevar(n));

	impl(n) = Argmax(mun);
	disp(sprintf('n=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g mu0=%g', ...
		n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
		mun(impl(n)), mat2str(extras.scale(:,n),3), extras.cov0(n), ...
		extras.noisevar(n), extras.mu0(n))); 
    end
end
