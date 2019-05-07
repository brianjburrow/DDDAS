%{
z=HuangBranin(x,y)
This function is similar to the six-hump camel back function from Branin(1972),
except it has the sign on one of the terms flipped.  This sign flipping appears
as a typo in Huang et al 2006, although it appears that paper actually used the
original Branin function for its computations.
%}
function z=HuangBranin(x,y)
	if (nargin == 1)
		% Allow passing in a vector instead of 2 arguments.
		[nrows,ncols]=size(x);
		if nrows == 2
			y = x(2,:);
			x = x(1,:);
		elseif ncols == 2
			y = x(:,2);
			x = x(:,1);
		else
			error('Single input argument of the wrong dimension'); 
		end
	end
	z = 4*x.^2 - 2.1*x.^4 + (1/3)*x.^6 + x.*y + 4*y.^2 + 4*y.^4;
	% Original Branin has the opposite sign on 4y^2. 
	%z = 4*x.^2 - 2.1*x.^4 + (1/3)*x.^6 + x.*y - 4*y.^2 + 4*y.^4;
end

