% [cov0hat, scalehat, betahat, dmodel, pv] = EstimateHyperparametersNoNoise(xd, y, scale0, scaleL, scaleU)
% Note this handles both isotropic & anisotropic cases.  Assumes noise variance
% is 0.
function [cov0hat, scalehat, betahat, dmodel, pv] = EstimateHyperparametersNoNoise(xd, y, scale0, scaleL, scaleU)
	n = length(y); % Number of measurements

	% dacefit rescales xd and y via (xd-mean(xd))./std(xd) and similarly
	% for y.  I had trouble figuring out whether they rescale all the
	% outputs, so I am going to scale the inputs myself, and then rescale
	% all the outputs.
	mxd = mean(xd);   sxd = std(xd);
	my = mean(y);   sy = std(y);
	xd = (xd - repmat(mxd,n,1)) ./ repmat(sxd,n,1);
	y = (y - my)/sy;

	% Also must rescale the initial guesses, and lower & upper bounds.
	% Scale is in units of xd, but is translation invariant, so it must be
	% rescaled only by std(xd).  The nugget g is unitless, and so does not
	% need to be rescaled.
	scale0 = scale0./sxd;
	scaleL = scaleL./sxd;
	scaleU = scaleU./sxd;

	% Transform coordinates to theta to be passed to corrgauss.
	theta0 = [scale0.^-2];
	thetaL = [scaleU.^-2];
	thetaU = [scaleL.^-2];

	% Retry multiple times to get a better solution.
	ntries = 10;
	theta=theta0;
	for i=1:ntries
		[dmodel,pv]=mydacefit(xd,y,@regpoly0,@corrgauss,theta,thetaL,thetaU);
		if (theta == dmodel.theta)
			break;
		end
		theta = dmodel.theta;
	end

	% Transform estimates back from dmodel
	cov0hat = dmodel.sigma2;
	k = length(scale0); % handles both isotropic, and anistropic cases.
	thetahat = dmodel.theta(1:k);
	scalehat = thetahat.^-.5;
	betahat = dmodel.beta;
	
	if (any(scalehat - scaleL < 1e-6))
		disp(sprintf('Estimate scalehat=%s and lower bound scaleL=%s are close.  The lower bound may need to be lowered.', mat2str(scalehat,3), mat2str(scaleL,3)));
	else if (any(scaleU - scalehat < 1e-6))
		disp(sprintf('Estimate scalehat=%s and upper bound scaleU=%s are close.  The upper bound may need to be raised.', mat2str(scalehat,3), mat2str(scaleU,3)));
	end

	% Rescale back from our original rescaling.  beta is in units of y, and
	% must be rescaled for translation.  cov0 is in units of y^2, and is
	% translation invariant.  scale is in units of x, and is translation
	% invariant.
	betahat = betahat*sy + my;
	cov0hat = cov0hat * sy^2;
	scalehat = scalehat .* sxd;
end
