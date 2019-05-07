%{
[impl, x, y, extras] = SimEGO(N,sample,dim,n0,p,MPerDim)

 Simulates the EGO algorithm through N samples, and collects at each sample
 time 1<=n<=N the measurement and implementation decision.  Adaptively
 estimates the hyperparameters after a first stage.  Also assumes 0 measurement
 noise.  Note that, although EGO chooses measurement decisions that are good
 for minimization, this function does a simple manipulation to make EGO do
 maximization instead.

 N -- number of measurements to take.
 sample(x) -- function that, given a (flattened) measurement location x, returns a
   noisy sample of the function value at that point.
 dim -- number of dimensions in the search domain.
 n0 -- number of measurements to put in the first stage.
 MPerDim -- number of discretizations in each dimension.  A scalar.

 impl -- vector of implementation decisions
 x -- vector of measurement decisions
 y -- vector of observations
 extras -- structure containing:
	extras.cov0 -- vector of estimated cov0.
	extras.scale -- array of estimated scales.
%}
function [impl,x,y,extras] = SimEGO(N,sample,dim,n0,p,MPerDim)
    if (p~=2)
	error('SimEGO is hard-coded to p=2');
    end
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras.cov0 = zeros(1,N);
    extras.scale = zeros(dim,N);
    toZd = @(z) ZToZd(z,MPerDim,dim);
    noisevar = 0;

    % Run a first stage using latin hypercube sampling
    mu0 = zeros(M,1);
    xd(1:n0,:) = DiscreteLHS(MPerDim,dim,n0);
    for n=1:n0
	x(n) = ZdToZ(xd(n,:),MPerDim);
	y(n) = sample(x(n));
	if (n<n0)
		impl(n) = x(argmax(y(1:n)));
		extras.scale(:,n) = NaN;
		extras.cov0(n) = NaN;
		disp(sprintf('n=%d: EGO x=%d (%s)', n, x(n), mat2str(xd(n,:))));
	end
    end

    [mun,Sigma,cov0,scale] = helperUpdate(x,xd,y,M,n,p,MPerDim);
    extras.cov0(n)=cov0;
    extras.scale(:,n)=scale;
    impl(n) = argmax(mun);
    disp(sprintf('n=n0=%d: EGO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g noisevar=%g', ...
			    n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
				mun(impl(n)), mat2str(scale,3), cov0, noisevar)); 



    % Run the rest of the measurements
    for n=[n0+1:N]
	% mun, Sigma have been set up by the previous iteration of this loop.

	assert(all(diag(Sigma)>=0));
	sposterior = sqrt(diag(Sigma));

	% Because of numerical imprecision, we may have that sposterior is
	% slightly bigger than 0, even at alternatives we have measured.  We
	% explicitly set them to 0 here, although it is a hack.  Note that we
	% zero them in BayesUpdateBlock, but there is also one measurement
	% updated that occurs in NoninformativeBayesUpdate, which doesn't zero
	% them.
	sposterior(x(1:n-1))=0; 
	assert(all(sposterior(x(1:n-1))==0));
    	x(n) = EGO(-mun,sposterior); % flip the sign of mu to get EGO to maximize.
	xd(n,:) = toZd(x(n));
	y(n) = sample(x(n));

	% Debugging
	%{
	lambda = noisevar * zeros(M,1);
    	xckg = CorrelatedNormalKGWithSigma(mun,diag(diag(Sigma)),lambda);
	if (xckg ~= x(n))
		disp(sprintf('CKG(mun,diag(diag(Sigma))) disagrees with EGO, xckg=%d xego=%d',  xckg, x(n)));
		disp(diag(Sigma)');
		disp(mun');
	end
	%}

	% Do the update & get the implementation decision.
	[mun,Sigma,cov0,scale] = helperUpdate(x,xd,y,M,n,p,MPerDim);
	extras.cov0(n)=cov0;
	extras.scale(:,n)=scale;
	impl(n) = argmax(mun);
	disp(sprintf('n=%d: EGO x=%d (%s) impl=%d (%s) mun[impl]=%f scale=%s cov0=%g', ...
			n, x(n), mat2str(xd(n,:)), impl(n), mat2str(toZd(impl(n))), ...
			mun(impl(n)), mat2str(scale,3), cov0)); 
    end
end

function [mun, Sigma, cov0, scale] = helperUpdate(x,xd,y,M,n,p,MPerDim)
    % Calculate the covariance Sigma of our current belief, as a function
    % of the cov0 and scale, which specify Sigma at time 0, and the points
    % we have measured thus far.
    [cov0, scale] = EstimateHyperparametersNoNoise(xd(1:n,:),y(1:n));
    noisevar = 0;
    A0=PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
    B0 = ones(M);
    mu0 = zeros(M,1);

    % Do the updates with a regularization term, pretending that this is
    % the noise of the observation noise in order to avoid problems with
    % inverting an ill-conditioned matrix, and then fix the mu and Sigma.
    regularization = .015;
    [mun,Sigma]=NoninformativeBayesUpdate(mu0,A0,B0,x(1),regularization,y(1));
    noises = regularization*ones(size(x(2:n)));
    [mun,Sigma]=BayesUpdateBlock(mun,Sigma,x(2:n)',y(2:n)',noises');
    mun(x(1:n))=y(1:n);
    Sigma(x(1:n),:)=0;
    Sigma(:,x(1:n))=0;
end
