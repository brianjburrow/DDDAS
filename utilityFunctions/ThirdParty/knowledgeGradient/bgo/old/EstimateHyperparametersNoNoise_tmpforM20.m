function [sigma2, scale, theta] = EstimateHyperparametersNoNoise(xd, y, corr, theta0, lob, upb)

	% dacefit will give an error if xd has duplicates within it.  Remove
	% these duplicates.  
	% PF: I should also check that the y values for the removed duplicates
	% are the same, but I do not.
	[uxd,i,j] = unique(xd,'rows');
    y = y(i);
    xd = uxd;
	
	if (nargin == 2)
		corr = @corrgauss;
		theta0 = 10;
		lob = 1e-4;
		upb = 1e3;
	end
	if (nargin == 3)
		theta0 = 10;
		lob = 1e-4;
		upb = 1e3;
	end

	% The inputted theta0,lob,upd, as well as the fitted theta returned by
	% dacefit are for a rescaled version of xd and y.  The rescaling they
	% use is (xd-mean(xd))./std(xd).  The estimated theta they return is
	% for the rescaled distances.  The estimated sigma2 they return,
	% however, is for the original data, so this does not need to be
	% rescaled. 

	theta0 = theta0 .* std(xd).^2;
	lob = lob .* std(xd).^2;
	upb = upb .* std(xd).^2;

	ntries = 10;
	for i=1:ntries
		dmodel=dacefit(xd,y,@regpoly0,corr,theta0,lob,upb);
		theta0 = dmodel.theta;
	end

	xstd = dmodel.Ssc(2,:); % Can be a vector in the anisotropic case.
	sigma2 = dmodel.sigma2;
	theta = dmodel.theta ./ (xstd.^2);
	scale = 1./theta;

	if (theta - lob < 10e-6)
		warning('estimated theta is very close to lower bound.  The lower bound may need to be lowered.');
	else if (upb - theta < 10e-6)
		warning('estimated theta is very close to upper bound.  The upper bound may need to be raised.');
	end
end
