% [newmu, newSigma] = BayesUpdateBlockZero(mu, Sigma, x, y)
% The inputs x and y are vectors.  The function assumes that all measurements
% are perfect.
function [newmu, newSigma] = BayesUpdateBlockZero(mu, Sigma, x, y)

% Preprocess x,y to get rid of measurements for which we already know the
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
	end
	x = uniq_x;
	y = uniq_y;
end
% Second part.
d = diag(Sigma);
known = find(d(x)==0);
if (length(known)>0)
	unknown = setdiff(1:length(x),known);
	x=x(unknown);
	y=y(unknown);
end
if (length(x)==0)
	newmu = mu;
	newSigma = Sigma;
	return;
end

% The preprocessing is done.  Here is the main body of the function.  We use
% the Kalman filter to do the necessary updating.  The notation (H,S,K) is from
% the Wikipedia entry on Kalman filtering.  When we calculate S we add on a
% regularization entry to keep the condition number from being too big.
% After the computation we then "fix" the result of the regularization by
% making the appropriate entries of newmu equal to the observations and the
% appropriate entries of newSigma 0.
regularization = .015;
perfect = find(diag(Sigma)==0);
perfectval = mu(perfect);
n = length(x); % number of measurements
M = length(mu); % size of state space
assert(all(size(Sigma)==[M M]));
assert(length(y)==n);
H = zeros(n,M);
H([1:n],x) = diag(ones(1,n)); % A fast way of doing "for i=1:n, H(i,x(i))=1".
S = H*Sigma*H' + regularization*eye(n);

% Fix ill-conditioning of S.
loop = 0;
added = regularization;
og_cond = cond(S);
while (cond(S)>1e4)
    cond_is_bad = 1;
    loop = loop+1;
    toadd = loop*regularization;
    S = S+toadd*eye(n);
    added = added + toadd;
end
if (loop > 0)
    disp(sprintf('Condition number was large, original cond(S)=%g, added regularization %g, new cond(S)=%g', og_cond, added, cond(S)));
end
assert(all(diag(S)>0)); % otherwise we will get NaN.

K = Sigma*H'*S^-1;
residual = y - H*mu;
newmu = mu + K*residual;
newmu(x) = y;
newmu(perfect) = perfectval;
if (nargout > 1)
	newSigma = (eye(M)-K*H)*Sigma;
	newSigma(x,:) = 0;
	newSigma(:,x) = 0;
	newSigma(perfect,:) = 0;
	newSigma(:,perfect) = 0;

	% Fix any diagonal entries that are negative
	mindiag = min(diag(newSigma));
	loop = 0;
	while (mindiag<0 && loop<10)
		newSigma = newSigma + (-mindiag+eps(-mindiag))*eye(M);
        	disp(sprintf('BayesUpdateBlock: fixing positive semi-definiteness of covariance matrix, adding %g, loop %d', mineig, loop));
		mindiag = min(diag(newSigma));
		loop = loop + 1;
    	end
    	if(mindiag<0)
		disp(sprintf('Even after fixing positive-semidefiniteness, Sigma has min(eig)=%g', mineig)); 
	end
	assert(all(diag(newSigma)>=0));
end
