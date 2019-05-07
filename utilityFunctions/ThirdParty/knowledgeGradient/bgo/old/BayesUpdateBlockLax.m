% [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises)
% The inputs x, y, and noises are vectors.
function [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises)

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
assert(all(diag(S)>0)); % otherwise we will get NaN.
K = Sigma*H'*S^-1;
residual = y - H*mu;
newmu = mu + K*residual;
if (nargout > 1)
	newSigma = (eye(M)-K*H)*Sigma;

	% This is a fix to handle numerical imprecision that occurs in the case
	% when lambda_x = 0.  Without this fix, we can sometimes get negative
	% values on the diagonal.  A better way might be to keep the square
	% root of the covariance matrix, as discussed in the literature.
	fix= find(diag(newSigma)<0);
	measured = x(find(noises==0));
	tozero = unique([fix; measured]);
	nfix = length(tozero)-length(measured);
	if (nfix > 0)
		% Then there are elements of fix not in measured.
		disp(sprintf('Setting Sigma(x,x)=0 for %d unmeasured alternatives',nfix));
	end
	k = length(tozero);
	if (k>0)
		newSigma(tozero,:) = zeros(k,M);
		newSigma(:,tozero) = zeros(M,k);
    	end

	% When the covariance matrix is badly scaled we sometimes have problems
	% with newmu not being equal to y at values we have measured.
	newmu(measured) = y(find(noises==0));

    %{
	mineig = min(eig(Sigma));
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
end
