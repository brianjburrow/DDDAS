% z=TiltedBranin(x,y)
% Returns the value of the tilted Branin function from Huang et al 2006,
% modified from Branin(1972).  
% 
% Huang et al. 2006 uses the domain x\in[-5,10], y\in[0,15] in which the
% function has 3 local minima, and 1 global minimum.  That paper claims the
% minimum value is -1.17, at x=-3.2, y=12.3, but it achieves -1.1855 at x=-3.2
% and y=12.4.
% 
% Test it via
% v = TiltedBranin(0,0);
% 	for x=[-5:.1:10]
%	  v = min(min(v,TiltedBranin(x,[0:.1:15])));
% 	end
% and then check that v is close to -1.17.
%
function z=TiltedBranin(x,y)
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
	a = 5.1/(4*pi^2);
	b = 5/pi;
	c = 10*(1-1/(8*pi));
	z = (y - a*x.^2 + b*x-6).^2 + c*cos(x) + 10 + 0.5*x;
end
