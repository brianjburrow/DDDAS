%[impl,x,y,extras] = SimIND(N,sample,dim,MPerDim,noisevar)
%
% Simulates Independent KG algorithm through N samples, and collects at each
% sample time 1<=n<=N the measurement and implementation decision.
%
% N -- number of measurements to take.
% sample(x) -- function that, given a (flattened) measurement location x, returns a
%   noisy sample of the function value at that point.
% dim -- number of dimensions in the search domain.
% MPerDim -- number of discretizations in each dimension.  A scalar.
%
% impl -- vector of implementation decisions
% x -- vector of measurement decisions
% y -- vector of observations
% extras -- NaN
%
function [impl,x,y,extras] = SimIND(N,sample,dim,MPerDim,noisevar)
    M = MPerDim^dim;
    impl = zeros(1,N); % preallocate for speed.
    x = zeros(1,N);
    y = zeros(1,N);
    extras = NaN;
    lambda = noisevar*ones(1,M);

    % We begin with a noninformative prior, and make it proper by measuring all
    % of {1,...M} once before measuring any twice.  Before the prior is proper,
    % the implementation choice argmax_x \mu_x isn't well defined because the
    % \mu_x of the unmeasured x are not well defined, as they are the means of
    % non-integrable distributions.  We define the implementation decision in
    % such cases to be the estimated best of those alternatives that have been
    % measured so far.
    mu = zeros(1,M);
    beta = zeros(1,M);
    for n=1:min(M,N)
	noninformative = find(beta == 0);
	x(n) = ChooseRandomElement(noninformative);
	y(n) = sample(x(n));
	mu(x(n)) = y(n);
	beta(x(n)) = 1/lambda(x(n));
	% To get the implementation decision, find the best among those
	% we have measured.
	impl(n) = x(Argmax(y(1:n)));
	%disp(sprintf('n=%d: IND x=%d', n, x(n)));
    end

    for n=M+1:N
    	x(n) = IndependentKG(mu',1./beta',lambda');
	y(n) = sample(x(n));
	noisebeta = 1/lambda(x(n));
	newbeta = noisebeta + beta(x(n));
	mu(x(n)) = (noisebeta*y(n) + beta(x(n))*mu(x(n)))/newbeta;
	beta(x(n)) = newbeta;
    	impl(n) = Argmax(mu);
	%disp(sprintf('n=%d: IND x=%d impl=%d mun[impl]=%f', n, x(n), impl(n), mu(impl(n))));
    end
    disp(sprintf('n=%d: IND x=%d impl=%d mun[impl]=%f', n, x(n), impl(n), mu(impl(n))));
end
