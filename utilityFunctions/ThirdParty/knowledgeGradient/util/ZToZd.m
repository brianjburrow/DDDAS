%zd = ZToZd(z, M, d)
%Generic discretization function.  Takes a point z between 1 and M^d, and
%returns a point in Z^d, with 0<=z_i<M.  Note that in Z, we have a minimum value
%of 1, so this is appropriate for indexing arrays in Matlab, but in Z^d our
%minimum values are 0, which is more natural for that space.
function zd = ZToZd(z, M, d)
	%{
	if (d==2)
		% Improve speed in the d=2 case.
		[zd(1),zd(2)] = ind2sub([M M], z);
		return;
	end
	%}
		
	% Remove the "artificial" 1 that was added to make it easy to index
	% arrays using z.
	z = z - 1; 

	if (round(d) ~= d || d <= 0);
		error('dimension must be a strictly positive integer');
	end
	if (z<0 || z>M^d-1)
		error('Input out of the domain');
	end
	for i=1:d
		zd(i) = mod(z,M);
		z = (z - zd(i)) / M;
	end
	assert(z == 0);
end

%{
Implementation note: originally I wanted to allow different amounts of
discretization in different dimensions.  I wanted to pass a vector M in Z^d,
where M_i is the maximum value for z_i.  The description was going to be,
"Takes a point M in Z^d, and an integer Z between 0 and \prod_i M_i, and
returns a point in Z^d with z_i between 0 and M_i."
%}
