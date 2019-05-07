%{
Plots the posterior belief and merit functions for the KG and EGO policies
that result from a sequence of (noise-free) measurements of a 1-dimensional
function.  Note that the plot we produce can be either for maximization or
minimization.

Inputs:
	trueVal: the underlying true function being measured.
	cov0,scale,mu0: hyperparameters of the prior (all scalars).
	x: vector of locations where we have already measured
	minimize: set to 1 if KG and EGO should minimize.  Set to 0 if they
		should maximize.  Default to 0.

Outputs:
	xckg: the point to measure under the KG policy
	xego: the point to measure under the EGO policy
		

 This function is based on a function that was used to create Figure 1 in
 	Frazier Powell and Dayanik (2009) The Knowledge-Gradient Policy for
 	Correlated Normal Beliefs. INFORMS Journal on Computing.  

%}
function [xckg, xego] = PlotEgoKg(trueVal,cov0,scale,mu0,x,minimize)

if (nargin<6)
	minimize = 0;
end

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
if (minimize)
	[xego,logei]=EGO(mun,stderr);
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(-mun,Sigma,zeros(size(mun)));
else
	[xego,logei]=EGO(-mun,stderr);
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(mun,Sigma,zeros(size(mun)));
end
clear Sigma

subplot(2,1,1)
plot(1:M,trueVal,'-k',x,y,'bo',1:M,mun,'-r',1:M,mun+2*stderr,'--r',1:M,mun-2*stderr,'--r')
xlabel('x');
ylabel('value');
legend('\theta', 'x^k, k<n', '\mu_x^n', '\mu^n_x+/-2\surd\Sigma^n_{xx}','Location','BestOutside');

subplot(2,1,2)
kg=exp(logkg);
ei=exp(logei);
plot(1:M,ei,'-g',1:M,kg,'-b',xego,ei(xego),'xg', xckg, kg(xckg),'xb')
xlabel('x');
ylabel('figure of merit');
legend('EGO', 'KG','Orientation','horizontal');
