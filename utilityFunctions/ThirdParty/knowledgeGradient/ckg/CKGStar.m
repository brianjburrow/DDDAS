function [xkg,maxlogslope,nstar,xkg1,maxlogqkg1]=CKGStar(mu, SigmaSqrt, noisevar, horizon)
M = length(mu);
%disp(sprintf('mu=%s', mat2str(mu,3)));
%disp(sprintf('diag(Sigma)=%s', mat2str(diag(Sigma),3)));
logslope = zeros(horizon,M);
for k=1:horizon % k is the number of effective measurements
  [tmp1,tmp2,logQ]=CorrelatedKG(mu,SigmaSqrt,noisevar/k);
  logslope(k,:)=logQ - log(k);
end

[maxlogslope_x, nstar_x] = max(logslope,[],1); % Best  each x.
xkg = Argmax(maxlogslope_x);
%plot(logslope(:,xkg)), pause
nstar = nstar_x(xkg);
maxlogslope = max(maxlogslope_x);
% Check for ties, if you wish.
tied = find(maxlogslope_x==maxlogslope);
ntied = length(tied);
if (ntied>1)
    %warning('CKGStar has a tie for alternative with the largest slope.  Choosing uniformly at random among them.');
    disp(sprintf('CKG has a tie.  best log slope = %g number tied=%d', maxlogslope, ntied));
    disp(sprintf('CKG tied alternatives: %s', mat2str(tied)));
    %disp(Sigma);
end

xkg1 = Argmax(logslope(1,:));
maxlogqkg1 = max(logslope(1,:));
if (maxlogslope_x(xkg1) < maxlogslope)
	disp(sprintf('CKG decision %d has max log slope %g, < best max log slope %g', xkg1, maxlogqkg1, maxlogslope));
end
