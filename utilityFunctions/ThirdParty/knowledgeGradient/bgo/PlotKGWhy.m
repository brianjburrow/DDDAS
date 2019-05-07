%{
Plots an "explanation" of why the KG policy wants to measure where it does.  In
particular, plots the posterior (subplot(3,1,1)), the value of the maximum of
the posterior as a function of the scaled observation z (subplot(3,1,2)), 
and the argmax of the posterior as a function of the scaled observation z
(subplot(3,1,3)).

Inputs:
	trueVal: the underlying true function being measured.
	cov0,scale,mu0: hyperparameters of the prior (all scalars).
	x: vector of locations where we have already measured
	x_examine: examine the KG factor for measurement at this location.
%}
function PlotKGWhy(trueVal,cov0,scale,mu0,x,x_examine)

M = length(trueVal);
y = trueVal(x);
[tmp,SigmaSqrt0]=PowExpCov(cov0,scale,2,M);
clear tmp;
mu0vec = mu0*ones(M,1);
[mun,SigmaSqrt]=BayesUpdateSqrtBlock(mu0vec,SigmaSqrt0,x,y,0);
clear SigmaSqrt0
Sigma = full(SigmaSqrt*SigmaSqrt');
clear SigmaSqrt
stderr=sqrt(diag(Sigma));
stderr(x)=0;
mun(x)=y;

[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(mun,Sigma,zeros(size(mun)));

subplot(3,1,3)
CKGMaxDist(mun,Sigma,zeros(1,M),x_examine)
