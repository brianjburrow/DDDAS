% Pass in x,y measured from a 1-dimensional function.  Also pass in
% hyperparameters cov0, scale, noisevar, and the true function.  This function
% will plot the true function with measurements, as well as the posterior mu
% and diag(Sigma).
function [mun,Sigma]=Test1DBayesUpdate(x,y,cov0,scale,noisevar,truth)
M=length(truth);
dim = 1;
assert(length(scale)==dim);
p = 2;

note=sprintf('cov0=%g scale=%g',cov0, scale);

mu0=zeros(M,1);
[Sigma0,SigmaSqrt0]=PowExpCov(cov0,scale,p,M);
[mun,SigmaSqrt]=BayesUpdateSqrtBlock(mu0,SigmaSqrt0,x,y,noisevar);
SigmaDiag = diag(SigmaSqrt*SigmaSqrt');
err = sqrt(SigmaDiag);

plot(1:M,truth,'g-',x,y,'bo',1:M,mun,'r-',1:M,mun+2*err,'r--',1:M,mun-2*err,'r--')
legend('truth','measurements','posterior mu','posterior mu +/- 2*stderr');
