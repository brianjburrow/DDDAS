% This function tests the runtime of CorrelatedKG.
function ok = TestTimeCKG()
	M=1000;
	lambda = ones(1,M);
	dim = 1;
	[mu0,Sigma0,SigmaSqrt0]=init1(M,1);
	tic, CorrelatedKG(mu0,SigmaSqrt0,lambda); runtime = toc;
	disp(sprintf('TestTimeCKG: test1(), runtime=%g seconds', runtime));
end

function [mu0,Sigma0,SigmaSqrt0]=test1(M,scale)
	mu0 = zeros(M,1);
	cov0 = 1;
	[Sigma0,SigmaSqrt0] = PowExpCov(cov0,scale,2,M,1);
	disp(sprintf('TestTimeCKG: test1() cov0=%f scale=%f M=%d mu0=zeros(M,1)', cov0, scale, M));
end

%{
function [mu0,Sigma0]=test2(MPerDim)
	dim = 2;
	M = MPerDim^dim;
	mu0 = 10*rand(M,1);
	cov0 = 10^(5*rand());
	scale = 10^(5*rand()-3);
	Sigma0 = PowerExponentialCovarianceMatrix(cov0,scale*ones(1,dim),dim,MPerDim);
end

function [mun,Sigma]=test3(MPerDim)
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
	disp(sprintf('TestTimeCKG: test3 n0=%d cov0=%f scale=%f dim=%d MPerDim=%d', n0, cov0, scale, dim, MPerDim));
	disp(sprintf('mu0=%s', mat2str(mu0,3)));
end
%}
