% Compute a power exponential covariance matrix using two different functions,
% PowExpCov, and an old simpler but deprecated function, and sees whether they
% agree.
function ok = TestPowExpCov()

cov0=5e6;
scale=[40 160];
noisevar=1;
MPerDim=21;
dim = 2;
p = 2;

[Sigma0,SigmaSqrt0]=PowExpCov(cov0,scale,p,MPerDim);
ok = 1;
err=max(max(abs(Sigma0-full(tril(SigmaSqrt0*SigmaSqrt0')))));
if (err>1e-6)
	disp(sprintf('TestPowExpCov: max(abs(Sigma0-SigmaSqrt0^2=%g))',err));
	ok = 0;
end

Sigma0_old = PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
err=max(max(abs(tril(Sigma0_old)-Sigma0)));
if (err>1e-6)
	disp(sprintf('TestPowExpCov: max(abs(Sigma0_old-Sigma0))=%g',err));
	ok = 0;
end


if (ok)
	disp('TestPowExpCov: OK');
else
	disp('TestPowExpCov: FAIL');
end
