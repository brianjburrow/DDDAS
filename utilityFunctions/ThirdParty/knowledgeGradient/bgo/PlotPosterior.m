%{
Plots the posterior belief on a 1-dimensional function.  This posterior belief
results from a sequence of measurements of a 1-dimensional function.  

Inputs:
	trueVal: the underlying true function being measured.
	cov0,scale,mu0: hyperparameters of the prior (all scalars).
	noisevar: measurement variance (scalar)
	x: vector of locations where we have already measured
	y: vector of values observed at the locations x. 
	belief_axis: argument to the "axis" command for the upper
		subplot on which the posterior belief is plotted.  Optional.

Outputs:
	xckg: the point to measure under the KG policy

%}
function xckg = PlotKg(trueVal,cov0,scale,mu0,noisevar,x,y,belief_axis)

if (nargin<8)
	belief_axis = 0;
end

M = length(trueVal);
[tmp,SigmaSqrt0]=PowExpCov(cov0,scale,2,M);
clear tmp;
mu0vec = mu0*ones(M,1);
[mun,SigmaSqrt]=BayesUpdateSqrtBlock(mu0vec,SigmaSqrt0,x,y,noisevar);
clear SigmaSqrt0
Sigma = full(SigmaSqrt*SigmaSqrt');
clear SigmaSqrt
stderr=sqrt(diag(Sigma));
minimize = 0;
if (minimize)
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(-mun,Sigma,zeros(size(mun)));
else
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(mun,Sigma,zeros(size(mun)));
end
clear Sigma

% [istar_val, istar] = max(mun);
% [true_istar_val, true_istar] = max(trueVal);
plot(1:M,trueVal,'-k',x,y,'bo',1:M,mun,'-r',1:M,mun+2*stderr,'--r',1:M,mun-2*stderr,'--r');

if (length(belief_axis)>1)
  axis(belief_axis);
end
xlabel('x');
ylabel('value');
% legend('\theta', 'x^k, k<n', '\mu_x^n', '\mu^n_x+/-2\surd\Sigma^n_{xx}','Location','BestOutside');
% legend('\theta', 'x^k, k<n', '\mu_x^n', '\mu^n_x+/-2\surd\Sigma^n_{xx}','Orientation','horizontal','Location','NorthOutside');
