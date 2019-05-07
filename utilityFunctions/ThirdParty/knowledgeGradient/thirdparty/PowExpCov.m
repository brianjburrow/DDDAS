% [R,C] = PowExpCov(cov0,scale,p,MPerDim,makefull)
% Returns the covariance matrix (R) and its Cholesky decomposition (C).  Both
% are returned as lower-diagonal sparse matrices.  Note that R is really
% symmetric, while C is really lower-diagonal.
%
% This implementation is more numerically stable and faster than the one I had
% originally.
%
% The implementation is heavily based upon code found in the package dacefit,
% which is copyrighted.  It also relies upon the function ind2sub2, which is a
% 3rd-party function.  Both of these are problematic for distribution.
%
% The "makefull" argument is optional, but if passed as true it will make the R
% returned as a full, symmetric, matrix.
function [R,C] = PowExpCov(cov0,scale,p,MPerDim,makefull)
assert(p==2);
theta = scale.^-2;
dim = length(scale);
m = MPerDim^dim;

% Calculate distances D between points
mzmax = m*(m-1) / 2;        % number of non-zero distances
ij = zeros(mzmax, 2);       % initialize matrix with indices
D = zeros(mzmax, dim);      % initialize matrix with distances
sz = repmat(MPerDim, 1, dim);
ll = 0;
for k = 1 : m-1
  ll = ll(end) + (1 : m-k);
  ij(ll,:) = [repmat(k, m-k, 1) (k+1 : m)']; % indices for sparse matrix
  iz = ind2sub2(sz,k);
  jz = ind2sub2(sz,ij(ll,2)); % ij(ll,2) is j.
  D(ll,:) = repmat(iz, m-k, 1) - jz; % differences between points

  % Also could have done, although it might have been less efficient,
  % since ij(ll,1) is i (and is identically k).
  % iz = ind2sub2(sz,ij(ll,1));
  % jz = ind2sub2(sz,ij(ll,2));
  % D(ll,:) = iz - jz;
end

r = corrgauss(theta, D);
idx = find(r > 0);   o = (1 : m)';   
mu = (10+m)*eps;
R = sparse([ij(idx,1); o], [ij(idx,2); o], [r(idx); ones(m,1)+mu]);
% Make R lower diagonal, instead of upper diagonal.
R = R';
% Cholesky factorization with check for pos. def.
if (nargout>1)
	[C rd] = chol(R,'lower');
	if  rd, error('not positive definite'); end
	C = sqrt(cov0)*C;
end
if (nargin>=5 && makefull)
	R=full(R+R'-diag(diag(R)));
end
R = cov0*R;
