% This compares the performance of IndKGStar_TangentPoint and IndKG.  It
% doesn't really produce output with which we can verify whether IndKGStar_* is
% doing the correct thing.
function [m,e]=TestIndKGStar(nruns)
if (nargin<1)
  nruns = 2;
end
if (nruns < 2)
  error('TestIndKGStar does not work with less than 2 runs');
end
mu0 = zeros(20,1);
mu0(1) = 2;
sigmasq0 = ones(20,1);
noisevar = 1;
N = 500;
d = zeros(nruns,N);
for run=1:nruns
    truth = normrnd(mu0,sigmasq0);
    mu_KG = mu0;
    sigmasq_KG = sigmasq0;
    mu_KGE = mu0;
    sigmasq_KGE = sigmasq0;


    for n=1:N
	i_KG = argmax(mu_KG);
	i_KGE = argmax(mu_KGE);
	d(run,n) = truth(i_KGE)-truth(i_KG);

	x_KG = IndKG(mu_KG,sigmasq_KG,noisevar);
	x_KGE = IndKGStar_TangentPoint(mu_KGE,sigmasq_KGE,noisevar,N-n+1);

	y_KG = truth(x_KG) + normrnd(0,noisevar);
	y_KGE = truth(x_KGE) + normrnd(0,noisevar);

	mu_KG(x_KG) = (mu_KG(x_KG)/sigmasq_KG(x_KG) + y_KG/noisevar)/(1/sigmasq_KG(x_KG) + 1/noisevar);
	sigmasq_KG(x_KG) = 1/(1/sigmasq_KG(x_KG) + 1/noisevar);
	
	mu_KGE(x_KGE) = (mu_KGE(x_KGE)/sigmasq_KGE(x_KGE) + y_KGE/noisevar)/(1/sigmasq_KGE(x_KGE) + 1/noisevar);
	sigmasq_KGE(x_KGE) = 1/(1/sigmasq_KGE(x_KGE) + 1/noisevar);
    end
    m = mean(d);
    e = std(d) / sqrt(nruns);
end
