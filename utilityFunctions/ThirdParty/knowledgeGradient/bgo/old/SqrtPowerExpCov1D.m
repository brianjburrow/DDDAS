% For some reason this implementation, which I took from dacefit, is more
% numerically stable than the one I have in SqrtPowerExpCov1DNaive.   I don't
% know why.  Maybe it's because the diagonal has 1+(10+m)*eps rather than just
% 1 on it...
function [C,R] = SqrtPowerExpCov1D(cov0,scale,p,m)
assert(p==2);
assert(length(scale)==1);
theta = scale^-2;

S = [1:m]';
n = 1;

% Calculate distances D between points
mzmax = m*(m-1) / 2;        % number of non-zero distances
ij = zeros(mzmax, 2);       % initialize matrix with indices
D = zeros(mzmax, n);        % initialize matrix with distances
ll = 0;
for k = 1 : m-1
  ll = ll(end) + (1 : m-k);
  ij(ll,:) = [repmat(k, m-k, 1) (k+1 : m)']; % indices for sparse matrix
  D(ll,:) = repmat(S(k,:), m-k, 1) - S(k+1:m,:); % differences between points
end

r = corrgauss(theta, D);
idx = find(r > 0);   o = (1 : m)';   
mu = (10+m)*eps;
R = sparse([ij(idx,1); o], [ij(idx,2); o], ...
  [r(idx); ones(m,1)+mu]);  
% Cholesky factorization with check for pos. def.
[C rd] = chol(R);
if  rd
	error('not positive definite');
end
C = sqrt(cov0)*C';
