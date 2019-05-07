function trueValList=MakeTrueValList(scale,M,N)
	if (nargin < 3)
		N=1000;
	end
	mu0 = zeros(M,1);
	cov0 = 0.5;
	p = 2;
	Sigma0 = PowExpCov(cov0,scale,p,M,1)
	trueValList = zeros(N,M);
	for n=1:1000
		trueValList(n,:)=mvnrnd(mu0,Sigma0);
	end
end
