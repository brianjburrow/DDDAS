function TestEstimateHyperparameters
	Test2();
end

function Test4a()
    logthetas=[-3:.1:3];
    for i=1:length(logthetas)
        theta = 10^logthetas(i);
        [dmodel_theta(i), est_theta(i)] = Test4(theta);
    end
    plot(logthetas,log10(dmodel_theta),logthetas,log10(est_theta),logthetas,logthetas);
    xlabel('log(true theta)');
    ylabel('log(estimated theta)');
    legend('dmodel, good theta0', 'EstimateHyperparametersNoNoise, default theta0', 'truth');
end

function Test3()
	scale = .2;
	cov0 = 2;
	N = 100;
	for i=1:100
		[estcov(i), estscale(i)] = Test1(cov0,scale);
	end
	subplot(2,1,1);
	hist(estscale);
	sort(estscale)'
	title('Histogram of scale estimates, true scale=.2, true cov0=2');
	subplot(2,1,2);
	hist(estcov);
	title('Histogram of cov0 estimates, true scale=.2, true cov0=2');
end

function Test2()
  %scales= [.01:.01:.15];
  scales= [.05:.05:2];
  for i=1:length(scales)
	[estimatedCov0(i), estimatedScale(i), trueCov0(i), trueScale(i)] = Test1(2,scales(i));
	disp(sprintf('cov0 estimate=%g truth=%g\t scale estimate=%s truth=%s', estimatedCov0(i), trueCov0(i), mat2str(estimatedScale(i)), mat2str(trueScale(i))));
  end
  subplot(2,1,1);
  plot(trueScale,log10(estimatedScale), trueScale,log10(trueScale));
  title('Performance of MLE, with true cov0=2, dim=1, and # points = 10');
  xlabel('true Scale');
  ylabel('log10(estimated Scale)');
  subplot(2,1,2);
  plot(trueScale,estimatedCov0);
  xlabel('true Scale');
  ylabel('estimated Cov0');
end

% Create a power-exponential(p=2) covariance matrix using known parameters, and
% then generate a multivariate normal random variable Y from it.  Use
% EstimateHyperparametersNoNoise to then estimate the parameters of the
% covariance matrix from Y, and see whether it's close.
function [estimatedCov0, estimatedScale, trueCov0, trueScale] = Test1(trueCov0, trueScale)
	% Tweak these
	MPerDim = 10;
	dim = 1;

	if (nargin == 0)
		trueCov0 = rand();
		trueScale = 20*rand(1,dim);
	end
	p = 2;
	M = MPerDim^dim;
	Sigma = PowExpCov(trueCov0, trueScale, p, MPerDim, 1);
	mu = zeros(1,M);
	Y = mvnrnd(mu,Sigma);
	for i=1:M
		x(i,:) = ZToZd(i,MPerDim,dim);
	end
		
	scale0 = .1;
	scaleU = trueScale * 10;
	scaleL = trueScale / 10;
	[estimatedCov0, estimatedScale] = EstimateHyperparametersNoNoise(x, Y, scale0, scaleL, scaleU);
end


% Like Test1, in that it creates a known covariance matrix, generates a random
% variable Y from it, and then estimates the covariance matrix from the data Y.
% It differs in that it uses a random design for x, rather than a grid.
function [estimatedCov0, estimatedScale, trueCov0, trueScale] = Test1a(trueCov0, trueScale)
	% Tweak these
	M = 10;

	p = 2;
	if (nargin == 0)
		trueCov0 = rand();
		trueScale = 20*rand();
	end
	x = rand(M,1); % Random design.
	Sigma = zeros(M,M);
	for i=1:M
		for j=1:M
			d = abs(x(i)-x(j));
			Sigma(i,j) = trueCov0 * exp(-(d/trueScale)^2);
		end
	end

	mu = zeros(1,M);
	Y = mvnrnd(mu,Sigma);
		
	scale0 = .1;
	scaleU = trueScale * 10;
	scaleL = trueScale / 100;
	[estimatedCov0, estimatedScale] = EstimateHyperparametersNoNoise(x, Y, scale0, scaleL, scaleU);
end

% This is the same as Test1, but implemented in a slightly different way.
function [estimatedCov0, estimatedScale, trueCov0, trueScale] = Test1b(trueCov0, trueScale)
	% Tweak these
	M = 10;

	p = 2;
	if (nargin == 0)
		trueCov0 = rand();
		trueScale = 20*rand();
	end
	x = [1:M]';
	Sigma = zeros(M,M);
	for i=1:M
		for j=1:M
			d = abs(x(i)-x(j));
			Sigma(i,j) = trueCov0 * exp(-(d/trueScale)^2);
		end
	end

	mu = zeros(1,M);
	Y = mvnrnd(mu,Sigma);
		
	scale0 = .1;
	scaleU = trueScale * 10;
	scaleL = trueScale / 100;
	[estimatedCov0, estimatedScale] = EstimateHyperparametersNoNoise(x, Y, scale0, scaleL, scaleU);
end



function [dmodel_theta, est_theta] = Test4(theta)
	lob = theta/10;
	upb = theta*10;
	theta0 = theta;
	M = 20;
	corr = @corrgauss;

	x = [1:M];
	x = x/std(x);
	d = zeros(M^2,1);
	for i=1:M
		for j=1:M
            idx = sub2ind([M,M],i,j);
			d(idx)=x(i)-x(j);
		end
	end
	r = corr(theta,d);
    r = reshape(r,M,M);
	mu = zeros(1,M);
	y = mvnrnd(mu,r);
	dmodel = dacefit(x',y',@regpoly0,corr,theta0,lob,upb);
    dmodel_theta = dmodel.theta;
    [est_cov0,est_scale,est_theta]=EstimateHyperparametersNoNoise(x',y');
    disp(sprintf('truth=%f dmodel=%f EstHyperparameters=%f',theta,dmodel.theta,est_theta));
end



% Create a power-exponential(p=2) covariance matrix using known parameters, and
% then generate a multivariate normal random variable Z from it, and then add
% independent noise.  Use EstimateHyperparameters to then estimate the
% parameters used, and see whether it's close.
function [estCov0, estNoise, estScale, estBeta, trueCov0, trueScale, trueNoise, trueBeta] = Test5(trueCov0, trueScale, trueNoise, trueBeta)
	% Tweak these
	MPerDim = 10;
	dim = 1;

	if (nargin == 0)
		trueCov0 = rand();
		trueScale = 20*rand(1,dim);
		trueNoise = rand();
		trueBeta = rand();
	end
	p = 2;
	M = MPerDim^dim;
	Sigma = PowExpCov(trueCov0, trueScale, p, MPerDim,1) + trueNoise*eye(M);
	mu = zeros(1,M);
	Y = mvnrnd(mu,Sigma);
	for i=1:M
		x(i,:) = ZToZd(i,MPerDim,dim);
	end
		
	scale0 = trueScale;
	scaleL = trueScale ./ 10;
	scaleU = trueScale .* 10;
	g0 = .5;
	gL = 1e-10;
	gU = 1;
	[estCov0, estNoise, estScale, estBeta] = EstimateHyperparameters(x, Y, scale0, scaleL, scaleU, g0, gL, gU);
end

function Test5a()
  nsamples = 100;
  trueCov0 = 1;
  trueScale = 1;
  trueNoise = 10;
  trueBeta = 1;
  for i=1:nsamples
	[estCov0(i), estNoise(i), estScale(i), estBeta(i)] = Test5(trueCov0, trueScale, trueNoise, trueBeta);
  end

  subplot(2,2,1)
  hist(estCov0);
  xlabel('estimated cov0');
  ylabel('frequency');
  title(sprintf('trueCov0=%g trueScale=%g trueNoise=%g trueBeta=%g', trueCov0, trueScale, trueNoise, trueBeta));

  subplot(2,2,2)
  hist(estNoise);
  xlabel('estimated noise');
  ylabel('frequency');

  subplot(2,2,3)
  hist(estScale);
  xlabel('estimated scale');
  ylabel('frequency');

  subplot(2,2,4)
  hist(estBeta);
  xlabel('estimated beta');
  ylabel('frequency');
end

% The same as Test5a, but without noise.
function Test5b()
  nsamples = 100;
  trueCov0 = 1;
  trueScale = 1;
  for i=1:nsamples
	[estCov0(i), estScale(i)] = Test1(trueCov0, trueScale);
  end

  subplot(2,1,1)
  hist(estCov0);
  xlabel('estimated cov0');
  ylabel('frequency');
  title(sprintf('trueCov0=%g trueScale=%g trueNoise=0', trueCov0, trueScale));

  subplot(2,1,2)
  hist(estScale);
  xlabel('estimated scale');
  ylabel('frequency');
end



%{
function Test5c()
  %scales= [.01:.01:.15];
  scales= [.05:.05:2];
  for i=1:length(scales)
	[estimatedCov0(i), estimatedNoise(i), estimatedScale(i), trueCov0(i), trueScale(i), trueNoise(i)] = Test5(2,scales(i),0.5);
	disp(sprintf('cov0 estimate=%g truth=%g\t scale estimate=%s truth=%s', estimatedCov0(i), trueCov0(i), mat2str(estimatedScale(i)), mat2str(trueScale(i))));
  end
  subplot(2,1,1);
  plot(trueScale,log10(estimatedScale), trueScale,log10(trueScale));
  title('Performance of MLE, with true cov0=2, dim=1, and # points = 10');
  xlabel('true Scale');
  ylabel('log10(estimated Scale)');
  subplot(2,1,2);
  plot(trueScale,estimatedCov0);
  xlabel('true Scale');
  ylabel('estimated Cov0');
end
%}


% Same as Test5b, but compares Test1a using random design against Test1 using a fixed design.
function Test5d()
  nsamples = 100;
  trueCov0 = 1;
  trueScale = 1;
  for i=1:nsamples
	[estCov0_grid(i), estScale_grid(i)] = Test1b(trueCov0, trueScale);
	[estCov0_rand(i), estScale_rand(i)] = Test1a(trueCov0, trueScale);
  end

  subplot(2,2,1)
  hist(estCov0_grid);
  xlabel('estimated cov0, grid');
  ylabel('frequency');
  title(sprintf('trueCov0=%g trueScale=%g trueNoise=0', trueCov0, trueScale));

  subplot(2,2,2)
  hist(estScale_grid);
  xlabel('estimated scale, grid');
  ylabel('frequency');


  subplot(2,2,3)
  hist(estCov0_rand);
  xlabel('estimated cov0, rand');
  ylabel('frequency');
  title(sprintf('trueCov0=%g trueScale=%g trueNoise=0', trueCov0, trueScale));

  subplot(2,2,4)
  hist(estScale_rand);
  xlabel('estimated scale, rand');
  ylabel('frequency');

end


