%[impl, x, y, extras] = SimNoisyCKG(N,sample,dim,n0,p,MPerDim)
%
% Simulates the Correlated KG algorithm through N samples, and collects at each sample
% time 1<=n<=N the measurement and implementation decision.  Adaptively
% estimates the hyperparameters after a first stage, including the measurement
% noise.  
%
% N -- number of measurements to take.
% sample(x) -- function that, given a (flattened) measurement location x, returns a
%   noisy sample of the function value at that point.
% dim -- number of dimensions in the search domain.
% n0 -- number of measurements to put in the first stage.
% MPerDim -- number of discretizations in each dimension.  A scalar.
%
% impl -- vector of implementation decisions
% x -- vector of measurement decisions
% y -- vector of observations
% extras -- structure containing:
%	extras.mu0 -- vector of estimated mu0.
%	extras.cov0 -- vector of estimated cov0.
%	extras.scale -- vector of estimated scales.
%	extras.noisevar -- vector of estimated noise variances.
%	extras.logkgfactor -- vector of maximum logkgfactor at each iteration.
%
function [impl,x,y,extras] = SimNoisyCKG(N,sample,dim,n0,p,MPerDim,design)
    if (p~=2)
	error('SimCKG is hard-coded to p=2');
    end
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras.logkgfactor = zeros(1,N);
    extras.cov0 = zeros(1,N);
    extras.mu0 = zeros(1,N);
    extras.noisevar = zeros(1,N);
    extras.scale = zeros(dim,N);
    toZd = @(z) ZToZd(z,MPerDim,dim);

    % Run a first stage using latin hypercube sampling.
    mu0 = zeros(M,1);
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
		impl(n) = x(argmax(y(1:n)));
		extras.scale(:,n) = NaN;
		extras.mu0(n) = NaN;
		extras.cov0(n) = NaN;
		extras.noisevar(n)=NaN;
		disp(sprintf('n=%d: CKG x=%d (%s)', n, x(n), mat2str(xd(n,:))));
	end
	extras.logkgfactor(n) = NaN; % We don't compute it in the first stage
    end

    [mun,Sigma,cov0,scale,noisevar,mu0] = helperUpdate(x(1:n),xd(1:n,:),y(1:n),M,p,MPerDim);
    extras.cov0(n)=mu0;
    extras.cov0(n)=cov0;
    extras.scale(:,n)=scale;
    extras.noisevar(n) = noisevar;
    impl(n) = argmax(mun);
    disp(sprintf('n=n0=%d: CKG x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g mu0=%g noisevar=%g', ...
		n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
		mun(impl(n)), mat2str(scale,3), cov0, mu0, noisevar)); 

    % Run the rest of the measurements
    for n=[n0+1:N]
	% mun, Sigma have been set up by the previous iteration of this loop.
	assert(all(diag(Sigma)>=0));

	lambda = noisevar * zeros(M,1);
    	[x(n),extras.logkgfactor(n)] = CorrelatedNormalKGWithSigma(mun,Sigma,lambda);
	xd(n,:) = toZd(x(n));
	y(n) = sample(x(n));

	% Do the update & get the implementation decision.
	[mun,Sigma,cov0,scale,noisevar,mu0] = helperUpdate(x(1:n),xd(1:n,:),y(1:n),M,p,MPerDim);
	extras.mu0(n)=mu0;
	extras.cov0(n)=cov0;
	extras.scale(:,n)=scale;
	extras.noisevar(n) = noisevar;
	impl(n) = argmax(mun);
	disp(sprintf('n=%d: CKG x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g mu0=%g noisevar=%g log(kgfactor)=%g', ...
			n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
			mun(impl(n)), mat2str(scale,3), cov0, mu0, noisevar, extras.logkgfactor(n))); 
    end
end


function [mun, Sigma, cov0, scale, noisevar, mu0] = helperUpdate(x,xd,y,M,p,MPerDim)
	[cov0, noisevar, scale, mu0] = EstimateHyperparameters(xd,y);

	% Calculate the covariance Sigma of our current belief, as a function
	% of the cov0 and scale, which specify Sigma at time 0, and the points
	% we have measured thus far.
	mun = mu0*ones(M,1);
	noises = noisevar*ones(size(x));
	Sigma=PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
	[mun,Sigma]=BayesUpdateBlock(mun,Sigma,x',y',noises');
end
