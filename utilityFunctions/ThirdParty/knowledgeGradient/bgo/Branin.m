%{
z=Branin(x,y)
Returns the value of the six-hump camel back function from Branin(1972).  

A word about the domain often used with this function:
Huang et al. 2006 uses the domain x\in[-1.6,2.4], y\in[-.8,1.2] in which the
function has 6 local minima, and 2 global minima.  A bigger domain containing
more features of interest is x\in[-2,3] and y\in[-1,1.5].  I wonder why Huang
chose the smaller domain?  Maybe because the bigger domain has a huge ramp
going up on the right, and it destroys the GP?

If you want to quickly plot this function, run the following matlab commands:
	[x,y]=meshgrid(-2:.05:3,-1:.05:1.5);
	z=Branin(x,y);
	surf(x,y,z);
%}
function z=Branin(x,y)
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
	z = 4*x.^2 - 2.1*x.^4 + (1/3)*x.^6 + x.*y - 4*y.^2 + 4*y.^4;
end
