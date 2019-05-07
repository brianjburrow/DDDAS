%{
[x,logEI] = EGO(Yhat, s)
Figures out the measurement decision used by EGO (Efficient Global
Optimization) from Jones, Schonlau, Welch 1998.  Assumes the inputs are
vectors, not multi-dimensional arrays.  The goal of the measurement decision
generated is MINIMIZATION.  Assumes that the measurements are perfect, i.e., no
noise.

Inputs:
Yhat -- Yhat(x) is the current Bayes estimate of Y(x).  That is, it is the mean
	of the posterior distribution of Y(x).
s --	s(x) is the standard deviation of the posterior distribution of Y(x).
	Should be exactly 0 for those alternatives that have been measured.
%}
function [x,logei] = EGO(Yhat, s)
	if any(size(Yhat)~=size(s))
		error('Yhat and s must have the same size');
    end
    assert(isreal(s));

	% We linearize the multi-dimensional arrays and work with them instead.
	sz = size(Yhat);
	maxind = prod(sz);
	Yhatlin = Yhat(1:maxind); 
	slin = s(1:maxind);
	
	measured = find(s==0);
	unmeasured = setdiff([1:maxind], measured);
	if (length(measured) == 0)
		disp('Calling EGO when nothing has been measured yet.  Returning the exploitation decision.');
		x = Argmax(Yhatlin);
		logei = -Inf*ones(size(Yhat));
        	return;
	end
	if (length(unmeasured) == 0)
		disp('Calling EGO when everything has been measured.  Returning a random decision.');
		logei = -Inf*ones(size(Yhat));
		%x = Argmax(Yhatlin); % exploitation decision
		x = Argmax(logei); % Random decision.
		return
        end;

	% f_min is the best value among those that have been measured.  
	f_min = min(Yhatlin(measured));
		
	% Then calculate the log of the EI (expected improvement).
	% We have the relation,
	% EI = diff.*normcdf(z) + slin(unmeasured).*normpdf(z);
	%    = slin(unmeasured).*(z.*normcdf(z) + normpdf(z));
	% but this introduces numerical instabilities, so instead we use the
	% LogEI function.
	diff = f_min - Yhatlin(unmeasured);
	z = diff./slin(unmeasured);
	logei(unmeasured) = log(slin(unmeasured)) + LogEI(z)';
	logei(measured) = -Inf;

	% Debugging
	%{
	logei
	subplot(2,1,1);
	plot(logei);
	subplot(2,1,2);
	plot(1:length(Yhatlin),Yhatlin-2*slin,'b',1:length(Yhatlin),Yhatlin+2*slin,'b');
	%}

	% Return the point with the largest EI.
	x = Argmax(logei); 
	if (length(find(logei==max(logei)))>1)
		%warning('EGO has a tie for alternative with the largest expected improvement.  Choosing uniformly at random among them.');
		disp(sprintf('EGO has a tie.  max(logEI)=%g length(find(logEI==max(logEI)))=%d', max(logei), length(find(logei==max(logei)))));
		disp(sprintf('EGO tied alternatives: %s', mat2str(find(logei==max(logei)))));
		%disp(Sigma);
	end

	ismeasured = length(find(measured==x));
	if (ismeasured)
		% PF: This can happen because we are not careful in controlling
		% numerical precision when we compute EI.  This also happens in
		% our implementations of SKO, and Correlated KG.  We should
		% really be computing the logarithm of the EI, and using
		% asymptotic approximations to Mill's ratio where appropriate. 
		disp('EGO is measuring a known alternative although unknown alternatives remain.'); 
	end
end
