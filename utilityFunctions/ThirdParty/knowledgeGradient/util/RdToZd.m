%{
z = RdToZd(x, delta, origin)
Generic discretization function.  Takes a point in R^d, a discretization size
delta, and the location of the origin in R^d, and returns a discretized point
in Z^d.  It rounds toward the nearest value, rather than to +infinity or
-infinity.  delta can either be a scalar, or a vector of the same length as x.
origin is an optional argument, and if not supplied will be assumed to be 0 (in
R^d).
%}
function z = RdToZd(x, delta, origin)
	if (nargin == 2)
		origin = zeros(size(x));
	end
	if (any(size(x) ~= size(origin)))
		error('x and origin should have the same size');
	end
	if (length(delta)~= 1 && any(size(x) ~= size(delta)))
		error('delta must either be a scalar, or the same size as x');
	end
	z = round((x - origin) ./ delta);
end
