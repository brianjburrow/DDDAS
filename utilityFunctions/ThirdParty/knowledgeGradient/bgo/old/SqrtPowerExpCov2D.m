function SqrtSigma = SqrtPowerExpCov2D(cov0, scale, p, MPerDim)
assert(length(scale)==2);
M = MPerDim^2;
%Sigma = spalloc(M,M,M*(M-1)); % M*(M-1) non-zeros is enough to fill the upper diagonal
Sigma = zeros(M,M);
for i=1:M
    j = 1:M;
    [i1zd i2zd] = ind2sub([MPerDim MPerDim],i); % i is a single 2-dim point
    [j1zd j2zd] = ind2sub([MPerDim MPerDim],j); % vector of 2-dim points j
    dist1 = abs(i1zd-j1zd)/scale(1); % distance in dimension 1 between i and j
    dist2 = abs(i2zd-j2zd)/scale(2); % distance in dimension 2 between i and j
    Sigma(j,i) = exp(-(dist1.^p + dist2.^p)); % Just fill in upper diagonal
end
% We use the incomplete cholesky decomposition because sometimes Sigma is not
% positive semi-definite.
%SqrtSigma = sqrt(cov0) * cholinc(sparse(Sigma),1e-15)';
SqrtSigma = sqrt(cov0) * chol(Sigma,'lower');
