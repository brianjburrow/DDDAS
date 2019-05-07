function [sigma2Zhat, sigma2Ehat, scalehat, betahat, dmodel, pv] = EstimateHyperparameters(xd, y, theta0, lob, upb, ntries)
	if (nargin == 2)
		theta0 = 10;
		lob = 1e-4;
		upb = 1e3;
		ntries = 1;
	end

	% The inputted theta0,lob,upd, as well as the fitted theta returned by
	% dacefit are for a rescaled version of xd and y.  The rescaling they
	% use is (xd-mean(xd))./std(xd) and similarly for y.  The estimated
	% theta and beta they return is for the rescaled distances.  The
	% estimated sigma2 they return, however, is for the original data, so
	% this does not need to be rescaled. 

	theta0 = theta0 .* std(xd).^2;
	lob = lob .* std(xd).^2;
	upb = upb .* std(xd).^2;

	for i=1:ntries
		theta0noise = [theta0,.5];
		lobnoise = [lob,1e-10]; % Really, we would like 0 instead of 1e-10, but this is not allowed by dacefit.
		upbnoise = [upb,1];
		[dmodel,pv]=mydacefit(xd,y,@regpoly0,@corrgaussnoise,theta0noise,lobnoise,upbnoise);
		theta0noise = dmodel.theta;
	end

	xstd = dmodel.Ssc(2,:); % Can be a vector in the anisotropic case.
	sigma2hat = dmodel.sigma2;
	k = length(theta0); % Do this to handle both isotropic, and anistropic cases.
	ghat = dmodel.theta(k+1); % PF: does this need to be rescaled????
	thetahat = dmodel.theta(1:k) ./ (xstd.^2);
	scalehat = 1./thetahat;
	betahat = dmodel.beta*std(y) + mean(y);

	% sigma2hat = sigma2Zhat + sigma2Ehat,
	% and ghat = sigma2Zhat / (sigma2Zhat + sigma2Ehat), 
	% so we solve this and get the following:
	sigma2Zhat = ghat*sigma2hat;
	sigma2Ehat = sigma2hat - sigma2Zhat;

	if (thetahat - lob < 10e-6)
		warning('estimated theta is very close to lower bound.  The lower bound may need to be lowered.');
	else if (upb - thetahat < 10e-6)
		warning('estimated theta is very close to upper bound.  The upper bound may need to be raised.');
	end
end

