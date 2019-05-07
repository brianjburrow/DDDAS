%[impl, x, y, extras] = SimNoisySKO(N,sample,dim,n0,p,MPerDim)
%
% Simulates the SKO (Sequential Kriging Optimization) algorithm through N
% samples, and collects at each sample time 1<=n<=N the measurement and
% implementation decision.  Adaptively estimates the hyperparameters after a
% first stage, including the measurement noise.  
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
%	extras.cov0 -- vector of estimated cov0.
%	extras.scale -- array of estimated scales.
%	extras.noisevar -- array of estimated noise variances.
%
function [impl,x,y,extras] = SimNoisySKO(N,sample,dim,n0,p,MPerDim)
    if (p~=2)
	error('SimSKO is hard-coded to p=2');
    end
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras.cov0 = zeros(1,N);
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
		extras.cov0(n)=NaN;
		extras.noisevar(n)=NaN;
		extras.scale(:,n)=NaN;
		disp(sprintf('n=%d: SKO x=%d (%s)', n, x(n), mat2str(xd(n,:))));
	end
    end

    [mun,Sigma,cov0,scale,noisevar] = helperUpdate(x,xd,y,M,n,p,MPerDim);
    extras.cov0(n)=cov0;
    extras.scale(:,n)=scale;
    extras.noisevar(n) = noisevar;
    impl(n) = argmax(mun);
    disp(sprintf('n=n0=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g', ...
		n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
		mun(impl(n)), mat2str(scale,3), cov0, noisevar)); 

    % Run the rest of the measurements
    for n=[n0+1:N]
	% mun, Sigma have been set up by the previous iteration of this loop.
	assert(all(diag(Sigma)>=0));

	% Debugging check.  
	if (noisevar == 0)
	for i=2:n-1
		deltamu = abs(mun(x(i))-y(i));
		if (deltamu > 1e-10)
			disp(sprintf('abs(mun(x(i)-y(i)))=%g, i=%d', deltamu, i));
		end
		deltaSig = Sigma(x(i),x(i));
		if (deltaSig > 1e-10 || deltaSig < 0)
			disp(sprintf('Sigma(x(i),x(i))=%g, i=%d', deltaSig, i));
		end
		%assert(abs(mun(x(i))-y(i))<1e-9);
		%assert(Sigma(x(i),x(i))>=0 && Sigma(x(i),x(i))<1e-9);
	end
	end

	assert(all(diag(Sigma)>=0));
	sposterior = sqrt(diag(Sigma));
	c = 1; % optional argument to SKO
	x(n) = SKO(-mun,sposterior,sqrt(noisevar),c);
	xd(n,:) = toZd(x(n));
	y(n) = sample(x(n));

	% Do the update & get the implementation decision.
	[mun,Sigma,cov0,scale,noisevar] = helperUpdate(x,xd,y,M,n,p,MPerDim);
	extras.cov0(n)=cov0;
	extras.noisevar(n) = noisevar;
	extras.scale(:,n)=scale;
	impl(n) = argmax(mun);
	disp(sprintf('n=%d: SKO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g', ...
			n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
			mun(impl(n)), mat2str(scale,3), cov0, noisevar)); 
    end
end


function [mun, Sigma, cov0, scale, noisevar] = helperUpdate(x,xd,y,M,n,p,MPerDim)
	[cov0, noisevar, scale] = EstimateHyperparameters(xd(1:n,:),y(1:n));

	% Calculate the covariance Sigma of our current belief, as a function
	% of the cov0 and scale, which specify Sigma at time 0, and the points
	% we have measured thus far.
	A0=PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
	B0 = ones(M);
	mu0 = zeros(M,1);
	[mun,Sigma]=NoninformativeBayesUpdate(mu0,A0,B0,x(1),noisevar,y(1));
	noises = noisevar*ones(size(x(2:n)));
	[mun,Sigma]=BayesUpdateBlock(mun,Sigma,x(2:n)',y(2:n)',noises');
end
