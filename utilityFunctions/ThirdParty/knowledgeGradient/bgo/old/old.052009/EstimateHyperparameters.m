% [cov0hat, noisevarhat, scalehat, betahat, dmodel, pv] = EstimateHyperparameters(xd, y, scale0, scaleL, scaleU, g0, gL, gU)
% Note this handles both isotropic & anisotropic cases.  The nugget g is the
% ratio of process variance divided by (process variance + measurement noise
% variance), and so naturally has a lower bound of 0 and an upper bound of 1,
% but dacefit requires that its lower bound be strictly bigger than 0.  Just
% pick a small number like 1e-10.
function [cov0hat, noisevarhat, scalehat, betahat, dmodel, pv] = EstimateHyperparameters(xd, y, scale0, scaleL, scaleU, g0, gL, gU)
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

	% Transform coordinates to theta to be passed to corrgaussnoise_v2.
	% We work with -log1p(-g) instead of g because g is typically a number
	% like 1-10^-k, with k being a small integer.  When dacefit optimizes,
	% it works better on the scale of k than it does on the scale of g.  We
	% use the negative signs to make sure that the transformed g is
	% positive, which dacefit requires.  Note that if gU is 1 then
	% -log1p(-gU) is infinity.  This is ok.
	theta0 = [scale0.^-2 -log1p(-g0)];
	thetaL = [scaleU.^-2 -log1p(-gL)];
	thetaU = [scaleL.^-2 -log1p(-gU)];

	% Retry multiple times to get a better solution.
	ntries = 10;
	theta=theta0;
	for i=1:ntries
		[dmodel,pv]=mydacefit(xd,y,@regpoly0,@corrgaussnoise_v2,theta,thetaL,thetaU);
		if (theta == dmodel.theta)
			break;
		end
		theta = dmodel.theta;
	end

	% Transform estimates back from dmodel
	sigma2hat = dmodel.sigma2;
	k = length(scale0); % handles both isotropic, and anistropic cases.
	ghat = -expm1(-dmodel.theta(k+1));
	thetahat = dmodel.theta(1:k);
	scalehat = thetahat.^-.5;
	betahat = dmodel.beta;
	% sigma2hat = sigma2Zhat + sigma2Ehat,
	% and ghat = sigma2Zhat / (sigma2Zhat + sigma2Ehat), 
	% so we solve this and get the following:
	cov0hat = ghat*sigma2hat;
	noisevarhat = sigma2hat - cov0hat;
	
	if (any(scalehat - scaleL < 1e-6))
		disp(sprintf('Estimate scalehat=%s and lower bound scaleL=%s are close.  The lower bound may need to be lowered.', mat2str(scalehat,3), mat2str(scaleL,3)));
	else if (any(scaleU - scalehat < 1e-6))
		disp(sprintf('Estimate scalehat=%s and upper bound scaleU=%s are close.  The upper bound may need to be raised.', mat2str(scalehat,3), mat2str(scaleU,3)));
	end

	% Rescale back from our original rescaling.  beta is in units of y, and
	% must be rescaled for translation.  cov0 and noisevar are in units of
	% y^2, and are translation invariant.  scale is in units of x, and is
	% translation invariant.
	betahat = betahat*sy + my;
	cov0hat = cov0hat * sy^2;
	noisevarhat = noisevarhat * sy^2;
	scalehat = scalehat .* sxd;
end
