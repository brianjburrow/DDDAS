%{
x = ZdToRd(z, delta, origin)
Generic discretization function.  Takes a point in Z^d, a discretization size
delta, and the location of the origin in R^d, and returns an "un-discretized"
point in R^d.  delta can either be a scalar, or a vector of the same length as
z.  origin is an optional argument, and if not supplied will be assumed to be 0
(in R^d).

%}
function x = ZdToRd(z, delta, origin)
	if (nargin == 2)
		origin = zeros(size(z));
	end
	if (any(size(z) ~= size(origin)))
		error('z and origin should have the same size');
	end
	if (length(delta)~= 1 && any(size(z) ~= size(delta)))
		error('delta must either be a scalar, or the same size as z');
	end

	x = z.*delta + origin;
end
