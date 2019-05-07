% Pass in x,y measured from the Branin function with MPerDim=21 and the
% standard domain.  Also pass in hyperparameters cov0, scale, noisevar.
% This function will plot the true function with measurements, as well as the
% posterior mu and diag(Sigma).
function [mun,Sigma]=TestBraninBayesUpdate(x,y,cov0,scale,noisevar)
v=load('BayesUpdateTestData');
BraninVal=v.BraninVal;
MPerDim=21;
dim = 2;
M = MPerDim^dim;
p = 2;

note=sprintf('cov0=%g scale=%s',cov0, mat2str(scale,3));

mu0=ones(MPerDim^dim,1);
[Sigma0,SigmaSqrt0]=PowExpCov(cov0,scale,p,MPerDim);
[mun,SigmaSqrt]=BayesUpdateSqrtBlock(mu0,SigmaSqrt0,x,y,noisevar);
SigmaDiag = diag(SigmaSqrt*SigmaSqrt');
Sigmadiag_matrix = zeros(21,21);
for i=1:21
  for j=1:21
	ind=ZdToZ([i-1 j-1],21);
	Sigmadiag_matrix(i,j)=SigmaDiag(ind);
	mun_matrix(i,j)=mun(ind);
  end
end

% Compute xd
for n=1:length(x)
	xd(n,:) = ZToZd(x(n),MPerDim,dim);
end

subplot(2,3,1)
contour(BraninVal)
hold on
plot(xd(:,1)+normrnd(zeros(size(xd(:,1))),.1),xd(:,2)+normrnd(zeros(size(xd(:,2))),.1),'o')
title('Branin, with measurements');
hold off
subplot(2,3,4)
surfc(BraninVal)
title('Branin')

subplot(2,3,2)
contour(mun_matrix)
title(sprintf('mu, %s', note))
hold on
plot(xd(:,1)+normrnd(zeros(size(xd(:,1))),.1),xd(:,2)+normrnd(zeros(size(xd(:,2))),.1),'o')
hold off
subplot(2,3,5)
surfc(mun_matrix)
title('mu')

subplot(2,3,3)
contour(Sigmadiag_matrix)
title('Sigmadiag')
hold on
plot(xd(:,1)+normrnd(zeros(size(xd(:,1))),.1),xd(:,2)+normrnd(zeros(size(xd(:,2))),.1),'o')
hold off
subplot(2,3,6)
surfc(Sigmadiag_matrix)
title('Sigmadiag')




