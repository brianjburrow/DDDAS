%{
[x,logei] = SKO(Yhat, s, snoise, measured, c)
Figures out the measurement decision used by SKO (Sequential Kriging
Optimization) from Huang et al. 2006.  Returns the measurement decision as a
linear index.  The goal of the measurement decision generated is MINIMIZATION.

Inputs:
Yhat -- Yhat(x) is the current Bayes estimate of Y(x).  That is, it is the mean
	of the posterior distribution of Y(x).
s --	s(x) is the standard deviation of the posterior distribution of Y(x).
Both Yhat and s can be passed as multidimensional arrays.
snoise -- scalar giving the standard deviation of the noise.
measured -- vector of those indices that have already been measured.
c --	scalar parameter indicating risk aversion, for use in finding the
	"effective best point".  Optional argument; by default it is 1.
%}
function [x,logei] = SKO(Yhat, s, snoise, measured, c)
	if (snoise == 0)
		disp('SKO has been called with snoise==0.  Returning EGO decision.');
		[x,logei] = EGO(Yhat,s);
		return;
	end
		
	if (nargin<5)
		c = 1;
	end
	if any(size(Yhat)~=size(s))
		error('Yhat and s must have the same size');
	end

	% We linearize the multi-dimensional arrays and work with them instead.
	sz = size(Yhat);
	maxind = prod(sz);
	Yhatlin = Yhat(1:maxind); 
	slin = s(1:maxind);
	
	
	% If we haven't measured anything then xstar2 is not defined (see
	% below).  Then the whole algorithm is undefined.  We just return the
	% implementation decision.
	if (length(measured) == 0)
		warning('Calling SKO when nothing has been measured yet.  Returning the exploitation decision.');
		x = argmax(Yhatlin);
		% Neither EI nor Emax is not really defined.
		logei = -Info*ones(size(Yhat));
		return
	end


	% Find x**, which is the current "effective best point".  We linearize
	% u before taking the max, and keep x** as a linear index, rather than
	% a multi-dimensional subscript.  We are minimizing Yhat, so the
	% utility has negatives in it here, and then we maximize utility.
	u = -Yhat - c*s; 
	% I used to have the following bug, where I maximized u over the whole
	% domain.  SKO actually says to only take the maximum over those points
	% in the domain we have already measured.  The bug actually improved
	% performance a little bit.
	%[bestu, xstar2] = max(u(1:maxind));
	[bestu, xstar2] = max(u(measured)); % Note that we are maximizing utility.
	xstar2 = measured(xstar2); % Map it back to the correct indexing.
		
	% Then calculate EI (expected improvement).
	% We have the relation,
	% 	Emax = diff.*normcdf(z) + s.*normpdf(z);
	%	EI = Emax .* (1 - snoise./sqrt(s.^2 + snoise^2));
	% but this introduces numerical instabilities, so instead we use the
	% LogEI function.
	diff = Yhat(xstar2) - Yhat;
	z = diff./s;
	logemax = log(s) + LogEI(z)';
	logei = logemax + log1p(-snoise./sqrt(s.^2 + snoise^2));

	% Debugging
	%{
	subplot(3,1,1)
	plot(z,'o')
	subplot(3,1,2)
	Emax = diff.*normcdf(z) + s.*normpdf(z);
	EI = Emax .* (1 - snoise./sqrt(s.^2 + snoise^2));
	plot(1:length(EI),EI,1:length(EI),exp(logei))
	subplot(3,1,3)
	plot(1:length(Yhat),Yhat-2*s,'b',1:length(Yhat),Yhat+2*s,'b');
	%}
	
	% Return the point with the largest EI.
	[bestlogei, x] = max(logei(1:maxind)); 
end
