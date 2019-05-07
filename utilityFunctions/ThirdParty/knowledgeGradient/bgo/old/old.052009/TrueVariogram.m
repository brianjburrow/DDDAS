% Uses the exponential gaussian correlation function.
function vario=TrueVariogram(cov0,scale,M)
d = [0:M-1];
for i=1:length(d)
	% Originally I wrote it like this to make it more obvious how the
	% variogram is computed from the covariance.  
	% variogram(x1-x2) = Var(x1) + Var(x2) - 2*Cov(x1,x2)
	%vario(i)=(2*cov0-2*cov0*exp(-scale*d(i)^2));
	vario(i)=cov0*(1-exp(-scale*d(i)^2));
end
