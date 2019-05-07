%{
Like PlotEgoKg, but plots the posterior belief and merit function for just the
KG policy, and allows noise.  This posterior belief and merit function result
from a sequence of measurements of a 1-dimensional function.  Note
that the plot we produce can be either for maximization or minimization.

Inputs:
	trueVal: the underlying true function being measured.
	cov0,scale,mu0: hyperparameters of the prior (all scalars).
	noisevar: measurement variance (scalar)
	x: vector of locations where we have already measured
	y: vector of values observed at the locations x. 
	minimize: set to 1 if KG and EGO should minimize.  Set to 0 if they
		should maximize.  Default to 0.
	belief_axis: argument to the "axis" command for the upper
		subplot on which the posterior belief is plotted.  Optional.
	kg_axis: argument to the "axis" command for the lower subplot on which
		the kg factor is plotted.  Optional.

Outputs:
	xckg: the point to measure under the KG policy
		

 This function is based on a function that was used to create Figure 1 in
 	Frazier Powell and Dayanik (2009) The Knowledge-Gradient Policy for
 	Correlated Normal Beliefs. INFORMS Journal on Computing.  

%}
function xckg = PlotKg(trueVal,cov0,scale,mu0,noisevar,x,y,minimize,belief_axis,kg_axis)

if (nargin<8)
	minimize = 0;
	belief_axis = 0;
	kg_axis = 0;
end
if (nargin<9)
	belief_axis = 0;
	kg_axis = 0;
end
if (nargin<10)
	kg_axis = 0;
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
if (minimize)
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(-mun,Sigma,zeros(size(mun)));
else
	[xckg,maxlogkg,logkg]=CorrelatedNormalKGWithSigma(mun,Sigma,zeros(size(mun)));
end
clear Sigma

subplot(2,1,1)
[istar_val, istar] = max(mun);
[true_istar_val, true_istar] = max(trueVal);
plot(1:M,trueVal,'-k',x,y,'bo',1:M,mun,'-r',1:M,mun+2*stderr,'--r',1:M,mun-2*stderr,'--r', istar, istar_val, 'sr', true_istar, true_istar_val, 'sk');

if (length(belief_axis)>1)
  axis(belief_axis);
end
xlabel('x');
ylabel('value');
% legend('\theta', 'x^k, k<n', '\mu_x^n', '\mu^n_x+/-2\surd\Sigma^n_{xx}','Location','BestOutside');
% legend('\theta', 'x^k, k<n', '\mu_x^n', '\mu^n_x+/-2\surd\Sigma^n_{xx}','Orientation','horizontal','Location','NorthOutside');

subplot(2,1,2)
plot(1:M,logkg,'-b', xckg, logkg(xckg),'xb')
if (length(kg_axis)>1)
  axis(kg_axis);
end
xlabel('x');
ylabel('log(KG factor)');
% legend('KG factor', 'Location','BestOutside');
