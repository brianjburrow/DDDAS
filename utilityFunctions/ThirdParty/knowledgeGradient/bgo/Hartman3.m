% z=Hartman3(x,y)
% Returns the value of the Hartman 3 function from Hartman (1973).
% 
% Huang et al. 2006 uses the domain x,y,z\in[0,1], in which the function has
% many local minima, and 1 global minimum at x=.114,y=.556,z=.852 with a value
% of -3.86.
% 
% Test it via
% v = Hartman3(0,0,0);
% 	for x=[-5:.1:10]
%	  v = min(min(v,Hartman3(x,[0:.1:15])));
% 	end
% and then check that v is close to -3.86.
%
function ret=Hartman3(x,y,z)
	if (nargin == 1)
		% Allow passing in a vector instead of 3 arguments.
		[nrows,ncols]=size(x);
		if nrows == 3
			z = x(3,:);
			y = x(2,:);
			x = x(1,:);
		elseif ncols == 3
			z = x(:,3);
			y = x(:,2);
			x = x(:,1);
		else
			error('Single input argument of the wrong dimension'); 
		end
	end

	alpha = [3 10 30; .1 10 35; 3 10 30; .1 10 35];
	c = [1; 1.2; 3; 3.2];
	p = [.3689 .1170 .2673; .4699 .4387 .7470; .1091 .8732 .5547; .03815 .5743 .8828];

	ret = 0;
	for i=1:4
	  term =alpha(i,1)*(x-p(i,1)).^2 + ...
	  	alpha(i,2)*(y-p(i,2)).^2 + ...
		alpha(i,3)*(z-p(i,3)).^2;
	  ret = ret - c(i)*exp(-term);
	end
		
	% Here are other ways to do it, but I don't think it works with vectors x,y,z
	% term =	alpha(:,1)*(x-p(:,1)).^2 + ...
	% 	alpha(:,2)*(y-p(:,2)).^2 + ...
	% 	alpha(:,3)*(z-p(:,3)).^2;
	% term = alpha*((x-p).^2)';
	% res = c'*exp(-term);
end
