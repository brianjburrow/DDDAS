%z = ZdToZ(zd, M)
% Generic discretization function.  Takes a point z in Z^d and an integer M,
% where 0 <= z_i < M, and returns an integer betweer 1 and M^d.  Since the
% response begins at 1, it is appropriate for indexing into a matlab array.
function z = ZdToZ(zd, M)
	d = length(zd);
	if any(zd<0) || any(zd>M-1)
		error('The input zd violates the 0<=zd(i)<M constraint for some i.');
	end
	z = zd*(M.^[0:d-1])';
	assert(length(z)==1);
	z = z+1; % Make z begin at 1, appropriate for indexing a matlab array.
end



%{
Here is an old implementation that allowed M to be a vector, with M_i the
discretization level for dimension i, different for different dimensions.
function z = ZdToZ(zd, M)
	% Check that each z_i is between 0 and M_i.
	if any(zd<0) or any(zd>M)
		warn('Input out of the domain');
	end
	z = zd(1);
	for i=[2:length(M)]
		z.*M(i); % Like a register shift
		z = z + zd(i);
	end
end
%}
