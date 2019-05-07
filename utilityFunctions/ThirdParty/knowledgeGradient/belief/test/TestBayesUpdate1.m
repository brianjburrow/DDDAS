function ok = TestBayesUpdate1()

ok = 1;

cov0=5e6;
scale=[40 160];
load('BayesUpdateTestData.mat'); % Contains x, y, and BraninVal
noisevar=1;
MPerDim=21;
dim = 2;
p = 2;
N=67;

note=sprintf('cov0=%g scale=%s',cov0, mat2str(scale,3));

mu0=zeros(MPerDim^dim,1);
[Sigma0,SigmaSqrt0]=PowExpCov(cov0,scale,p,MPerDim);

% Recursive.  We don't do this anymore, since Batch is secretly recursive
% underneath, and doing it this way that is commented out is really slow
% because of all the cholesky decompositions.
%{
[mun,Sigma]=BayesUpdateBlock(mu0,Sigma0,x(1),y(1),1);
for n=2:67
	[mun,Sigma]=BayesUpdateBlock(mun,Sigma,x(n),y(n),noisevar);
end
for i=1:21
  for j=1:21
	ind=ZdToZ([i-1 j-1],MPerDim);
	mun_matrix(i,j)=mun(ind);
	Sigmadiag_matrix(i,j)=Sigma(ind,ind);
  end
end
%}

% Batch
[mun_batch,SigmaSqrt_batch]=BayesUpdateSqrtBlock(mu0,SigmaSqrt0,x',y',noisevar*ones(size(x')));
SigmaDiag_batch=diag(SigmaSqrt_batch*SigmaSqrt_batch');
for i=1:21
  for j=1:21
	ind=ZdToZ([i-1 j-1],21);
	Sigmadiag_batch_matrix(i,j)=SigmaDiag_batch(ind);
	mun_batch_matrix(i,j)=mun_batch(ind);
  end
end

% Compute xd
for n=1:length(x)
	xd(n,:) = ZToZd(x(n),MPerDim,dim);
end

subplot(2,3,1)
contour(BraninVal)
hold on
plot(xd(:,1)+normrnd(zeros(100,1),.1),xd(:,2)+normrnd(zeros(100,1),.1),'o')
title('Branin, with measurements');
hold off
subplot(2,3,4)
surfc(BraninVal)
title('Branin')

subplot(2,3,2)
contour(mun_batch_matrix)
title(sprintf('batch mu, %s', note))
hold on
plot(xd(:,1)+normrnd(zeros(100,1),.1),xd(:,2)+normrnd(zeros(100,1),.1),'o')
hold off
subplot(2,3,5)
surfc(mun_batch_matrix)
title('batch mu')

subplot(2,3,3)
contour(Sigmadiag_batch_matrix)
title('batch Sigmadiag')
hold on
plot(xd(:,1)+normrnd(zeros(100,1),.1),xd(:,2)+normrnd(zeros(100,1),.1),'o')
hold off
subplot(2,3,6)
surfc(Sigmadiag_batch_matrix)
title('batch Sigmadiag')




%{
subplot(2,3,1)
surfc(BraninVal)
title('Branin')
subplot(2,3,4)
contour(BraninVal)
hold on
plot(xd(:,1)+normrnd(zeros(100,1),.1),xd(:,2)+normrnd(zeros(100,1),.1),'o')
title('Branin, with measurements');
hold off

subplot(2,3,2)
surfc(mun_matrix)
title(sprintf('recursive mu, %s', note))
subplot(2,3,5)
diff = mun_matrix - mun_batch_matrix;
surfc(diff)
title(sprintf('recursive Sigma diag - batch Sigma diag, max diff=%g', max(max(abs(diff)))));

subplot(2,3,3)
surfc(Sigmadiag_matrix)
title('recursive Sigmadiag')
subplot(2,3,6)
diff = Sigmadiag_matrix - Sigmadiag_batch_matrix;
surfc(diff)
title(sprintf('recursive Sigma diag - batch Sigma diag, max diff=%g', max(max(abs(diff)))));
%}

