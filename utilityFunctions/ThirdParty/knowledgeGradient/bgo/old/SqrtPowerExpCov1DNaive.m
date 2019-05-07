function SqrtSigma = SqrtPowerExpCov1DNaive(cov0, scale, p, M)
assert(length(scale)==1);
% I thought that I would only have to fill in the upper diagonal, but
% apparently cholinc requires the entire matrix to be full.
% Sigma = spalloc(M,M,M*(M+1)/2); % M*(M+1)/2 non-zeros is enough to fill the upper diagonal
Sigma = zeros(M,M);
for i=1:M
    j = 1:M; % Would be 1:i if we wanted upper triangular only.
    d = abs(i-j);
    Sigma(j,i) = exp(-(d/scale).^p); % upper triangular if j=1:i;
end
% We use the incomplete cholesky decomposition because sometimes Sigma is not
% positive semi-definite.
%assert (max(max(abs(cov0*Sigma - PowerExponentialCovarianceMatrix1D(cov0,scale,p,M))))<1e-10);
SqrtSigma = sqrt(cov0)*cholinc(sparse(Sigma),1e-15)'; % lower triangular
% I also tried to use this...
%Sigma(1:M,1:M)=1+(10+M)*eps;
%SqrtSigma = sqrt(cov0)*chol(sparse(Sigma));
