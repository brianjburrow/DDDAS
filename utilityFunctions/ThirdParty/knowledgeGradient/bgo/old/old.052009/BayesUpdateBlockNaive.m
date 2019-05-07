% [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises, strict)
% The inputs x, y, and noises are vectors.
% This is the strict version of the function, and so it checks that certain
% properties of the posterior mu and Sigma are within tolerance.
function [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises, strict)
if (nargin == 5)
	% Make it medium-strict by default, i.e. not 0=lax, not 2=super-strict.
	strict = 0;
end

% Preprocess x,y,noises to get rid of measurements for which we already know the
% result perfectly.  If we do not do this, the matrix S will be singular and
% the inversion of it will lead to NaN.  This preprocessing consists of three
% parts.  First, if there are duplicates in x, we combine them together into
% one measurement.  Second, if we measure an alternative with the new x that we
% already knew perfectly under Sigma, we drop it.  Third, if the new number of
% measurements is zero, we just return the old mu and Sigma.
% First part.
uniq_x = unique(x);
if (length(uniq_x) < length(x))
	uniq_y = zeros(size(uniq_x));
	uniq_noises = zeros(size(uniq_x));
	for k=1:length(uniq_x)
		idx = find(x==uniq_x(k));
		uniq_y(k) = mean(y(idx));
		uniq_noises(k) = 1/sum(noises(idx).^-1); % works if noises are 0.
	end
	x = uniq_x;
	y = uniq_y;
	noises = uniq_noises;
end
% Second part.
d = diag(Sigma);
known = find(d(x)==0);
if (length(known)>0)
	unknown = setdiff(1:length(x),known);
	x=x(unknown);
	y=y(unknown);
	noises=noises(unknown);
end
if (length(x)==0)
	newmu = mu;
	newSigma = Sigma;
	return;
end

% The preprocessing is done.  Here is the main body of the function.  We use
% the Kalman filter to do the necessary updating.  The notation (H,S,K) is from
% the Wikipedia entry on Kalman filtering.
n = length(x); % number of measurements
M = length(mu); % size of state space
assert(all(size(Sigma)==[M M]));
assert(length(y)==n);
assert(length(noises)==n);
H = zeros(n,M);
H([1:n],x) = diag(ones(1,n)); % A fast way of doing "for i=1:n, H(i,x(i))=1".
S = H*Sigma*H' + diag(noises);
loop = 0;
while (cond(S)>1e5)
    %S = S+eye(n)*exp(-7+n); % I used to have this.
    S = S+eye(n)*exp(-7+loop);
    loop = loop+1;
    disp(sprintf('Condition number is large, cond(S)=%g, loop=%d', cond(S), loop));
end
assert(all(diag(S)>0)); % otherwise we will get NaN.
K = Sigma*H'*S^-1;
residual = y - H*mu;
newmu = mu + K*residual;
if (nargout > 1)
	% Maybe for not strict I should put fixes here on eye(M)-K*H?
	newSigma = (eye(M)-K*H)*Sigma;
end

% Post-processing checks and fixes, depending on strictness.
perfect = find(noises==0); % These are the perfect measurements.

% Check or fix newmu on perfect measurements.
if (strict>1)
	% Check that the posterior mu agrees with the measurements for perfect measurements.
	err = newmu(x(perfect)) - y(perfect);
	maxerr = max(abs(err));
	if (maxerr>1e-10)
		error(sprintf('newmu(x(perfect))=%s y(perfect)=%s maxerr=%g', mat2str(newmu(x(perfect))), mat2str(y(perfect)), maxerr));
	end
else
	% When the covariance matrix is badly scaled we sometimes have problems
	% with newmu not being equal to y at values we have measured.
	newmu(x(perfect)) = y(perfect);
end


if (nargout > 1 && strict>1)
	% Require positive semi-definiteness.
	mineig = min(eig(newSigma));
	if (mineig<0)
		error(sprintf('Sigma is not positive semi-definite, mineig=%g', mineig));
	end

	% Should I also put a check on the diagonals being 0?
elseif (nargout > 1) % not strict, fixes on Sigma
	% Handle negative definiteness.
	% This fix has problems because eig can be imaginary, or can be
	% negative but very close to 0, and we still get mineig=0 even though
	% the matrix is not positive semi-definite.
	%{
	mineig = min(eig(newSigma));
	loop = 0;
	while (mineig<0 && loop<10)
		newSigma = newSigma + (-mineig+eps(-mineig))*eye(M);
        	%disp(sprintf('BayesUpdateBlock: fixing positive semi-definiteness of covariance matrix, adding %g, loop %d', mineig, loop));
		mineig = min(eig(newSigma));
		loop = loop + 1;
    	end
    	if(mineig<0)
		disp(sprintf('Even after fixing positive-semidefiniteness, Sigma has min(eig)=%g', mineig)); 
	end
	%}

	mindiag = min(diag(newSigma));
	loop = 0;
	while (mindiag<0 && loop<10)
		newSigma = newSigma + (-mindiag+eps(-mindiag))*eye(M);
        	%disp(sprintf('BayesUpdateBlock: fixing positive semi-definiteness of covariance matrix, adding %g, loop %d', mineig, loop));
		mindiag = min(diag(newSigma));
		loop = loop + 1;
    	end
    	if(mindiag<0)
		disp(sprintf('Even after fixing positive-semidefiniteness, Sigma has min(eig)=%g', mineig)); 
	end

	% Set to 0 those entries corresponding to a perfect measurement.
	k = length(perfect);
	if (k>0)
		newSigma(x(perfect),:) = zeros(k,M);
		newSigma(:,x(perfect)) = zeros(M,k);
    	end


	%{
	% This is a fix to handle numerical imprecision that occurs in the case
	% when lambda_x = 0.  Without this fix, we can sometimes get negative
	% values on the diagonal.  
	fix= find(diag(newSigma)<0);
	tozero = unique([fix; perfect]);
	nfix = length(tozero)-length(perfect);
	if (nfix > 0)
		% Then there are elements of fix not in perfect.
		disp(sprintf('Setting Sigma(x,x)=0 for %d unmeasured alternatives',nfix));
	end
	k = length(tozero);
	if (k>0)
		newSigma(tozero,:) = zeros(k,M);
		Should it be x(tozero)?
		newSigma(:,tozero) = zeros(M,k);
    	end
	%}

	assert(all(diag(newSigma)>=0));
end
