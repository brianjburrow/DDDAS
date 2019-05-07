cov0=.5;
M=80;
scale=.5*(M-1);
N=20;
mu0 = zeros(M,1);
[Sigma0,SqrtSigma0] = PowExpCov(cov0,scale,2,M);
Sigma0 = SqrtSigma0*SqrtSigma0';
truth = mvnrnd(mu0,Sigma0);
x = DiscreteLHS(M,1,N)+1; % The +1 transforms from Z to index space.
noisevar=.01;
y = truth(x) + mvnrnd(zeros(size(x)),noisevar*eye(N))

[mun,SqrtSigma]=BayesUpdateSqrtBlock(mu0,SqrtSigma0,x',y',noisevar);
err = sqrt(diag(SqrtSigma*SqrtSigma'));

plot(1:M,truth,'g-',x,y,'bo',1:M,mun,'r-',1:M,mun+2*err,'r--',1:M,mun-2*err,'r--')
legend('truth','measurements','posterior mu','posterior mu +/- 2*stderr');
