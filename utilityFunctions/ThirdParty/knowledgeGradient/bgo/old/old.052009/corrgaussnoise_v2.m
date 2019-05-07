function  [r, dr] = corrgauss(theta, d)
% CORRGAUSSNOISE_V2  Gaussian correlation function, with measurement noise,
% and log1p scaling of the nugget g. 
%
%             n
%   r_i = g * prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%             j=1
%
% The parameter g is given by -expm1(-theta(n+1)).  g is the ratio of process
% variance divided by (process variance + measurement noise variance).  Fixing
% it to 1 gives the case with no measurement noise.
%
% If length(theta) = 2, then the model is isotropic, so theta_1,..,theta_n are
% given by theta(1).  In this case, g is given by theta(2).
%
% Call:    r = corrgauss(theta, d)
%          [r, dr] = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% pfrazier@princeton.edu
% Last update May 30, 2008
% Based on corrgauss.m from the DACE package written by hbn@imm.dtu.dk, version
% from June 2, 2002.

[m n] = size(d);  % number of differences and dimension of data
if  length(theta) == 2
  theta_d = repmat(theta(1),1,n);
  g = -expm1(-theta(2));
elseif  length(theta) == n+1
  theta_d = theta(1:n);
  g = -expm1(-theta(n+1));
else
  error(sprintf('Length of theta must be 1+1 or %d+1',n))
end

if (g<0 || g>1)
  error('g, given by the last element of theta, must satisfy 0<=g<=1');
end

td = d.^2 .* repmat(-theta_d(:).',m,1);
r = g*exp(sum(td, 2));

if  nargout > 1
  error('not supported');
  %the next line was the code from corrgauss, but I'm not sure what to do here
  %dr = repmat(-2*theta_d(:).',m,1) .* d .* repmat(r,1,n);
end
