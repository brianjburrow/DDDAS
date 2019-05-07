% This function computes argmax_{n\ge1} V(n)/n, where V(n) is the value of measuring an alternative n times.
% Think of this as the point where a line from (0,0) tangent to the curve (n,V(n)) touches the curve.
% delta and beliefvar can be a vectors of the same size, in which case this
% function will return a vector with the value computed for each pair.
% noisevar can either be a scalar, or a vector the same size as beliefvar.
% nmax and nmin should be scalars, containing the upper and lower locations to
% truncate.  These arguments are optional, and their default values are
% nmax=Infinity, and nmin=0.
function nstar = TangentPoint(delta, beliefvar, noisevar, nmin, nmax)
	if (nargin < 5)
		nmax = Inf;
	end
	if (nargin < 4)
		nmin = 0;
	end
	assert(length(nmax)==1);
	assert(length(nmin)==1);
	assert(nmax>=nmin);
	assert(all(size(beliefvar)==size(delta)));
	assert(length(noisevar)==1 || all(size(beliefvar)==size(noisevar)));
	n0 = noisevar ./ beliefvar;
	% This is old code, and there may be a bug in it, particularly when noisevar<>1.
	% t = noisevar - 2 + n0.*delta.^2;
	% nstar = (t+sqrt(8*n0.*delta.^2+t.^2)).*n0/4;
	u = delta.^2./beliefvar;
	nstar = (-1+u+sqrt(1+6*u+u.^2)).*n0/4;

	% Up to here, nstar is argmax_{n\ge0} V(n)/n.  We can truncate at a
	% minimum point, like 1, or a maximum point, which would be interpreted
	% as the number of points left.
	if (nmin > 0)
		nstar(find(nstar<nmin))=nmin;
	end
	if (nmax < Inf)
		nstar(find(nstar>nmax))=nmax;
	end
end
