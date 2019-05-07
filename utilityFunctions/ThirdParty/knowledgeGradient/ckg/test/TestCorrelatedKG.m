% This test does not return ok or not ok.  Instead, it produces a plot 

function TestCorrelatedKG()
	MPerDim=10;
	dim = 2;
	M = MPerDim^dim;
	[mu0,Sigma0]=init1(MPerDim, 1);
	t = 10.^[-2:.1:2];
	xkg = zeros(1,length(t));
	logkg = zeros(1,length(t));
	lambda = zeros(1,M);
	for i=1:length(t)
		[xkg(i),logkg(i)]=CorrelatedNormalKGWithSigma(mu0,t(i)*Sigma0,lambda);
	end

	% In the init1 case, the prediction is that kg at t should be sqrt(t) * kg at 1.
	logkg1 = logkg(find(t==1))
	predicted_kg = sqrt(t)*exp(logkg1);
	predicted_logkg = .5*log(t)+logkg1;

	subplot(2,1,1);
	plot(t,exp(logkg),'o',t,predicted_kg);
	legend('Numerical','Theoretical');
	ylabel('max(\nu^{KG}_x)');
	xlabel('t');

	subplot(2,1,2);
	plot(log(t),logkg,'o',log(t),predicted_logkg);
	legend('Numerical','Theoretical');
	ylabel('log(max(\nu^{KG}_x))');
	xlabel('log(t)');
end

function [mu0,Sigma0]=init1(MPerDim,scale)
	dim = 2;
	M = MPerDim^dim;
	mu0 = zeros(M,1);
	cov0 = 1;
	Sigma0 = PowerExponentialCovarianceMatrix(cov0,scale*ones(1,dim),dim,MPerDim);
	disp(sprintf('cov0=%f scale=%f dim=%d MPerDim=%d mu0=zeros(M,1)', cov0, scale, dim, MPerDim));
end

function [mu0,Sigma0]=init2(MPerDim)
	dim = 2;
	M = MPerDim^dim;
	mu0 = 10*rand(M,1);
	cov0 = 10^(5*rand());
	scale = 10^(5*rand()-3);
	Sigma0 = PowerExponentialCovarianceMatrix(cov0,scale*ones(1,dim),dim,MPerDim);
	disp(sprintf('cov0=%f scale=%f dim=%d MPerDim=%d', cov0, scale, dim, MPerDim));
	disp(sprintf('mu0=%s', mat2str(mu0,3)));
end

function [mun,Sigma]=init3(MPerDim)
	dim = 2;
	M = MPerDim^dim;
	mu0 = 10*rand(M,1);
	cov0 = 10^(5*rand());
	scale = 10^(5*rand()-3);
	Sigma0 = PowerExponentialCovarianceMatrix(cov0,scale*ones(1,dim),dim,MPerDim);
	n0 = MPerDim;
	xd = round(MPerDim*(lhsdesign(n0,dim,'smooth','off')-.5/n0));
	mun = mu0;
	Sigma = Sigma0;
	for n=1:n0
		x(n) = ZdToZ(xd(n,:),MPerDim);
		sample_mu = mun(x(n));
		sample_var = Sigma(x(n),x(n));
		y(n) = normrnd(sample_mu,sample_var);
		[mun,Sigma] = BayesUpdateBlock(mu0,Sigma0,x(1:n)',y(1:n)',zeros(1,n)');
	end
	disp(sprintf('init3 n0=%d cov0=%f scale=%f dim=%d MPerDim=%d', n0, cov0, scale, dim, MPerDim));
	disp(sprintf('mu0=%s', mat2str(mu0,3)));
end
