% This implementation of the independent normal KG(*) policy uses the
% TangentPoint code to approximate the value of nstar to use.
%
% mu should be an M-vector,
% beliefvar should be an M vector, containing belief variances.
% noisevar should either be a scalar homogeneous noise variance, or an M-vector of the noise variances.
% nleft should be a scalar, containing how many measurements are left to make.  Can be Infinity.
function [xkg,maxLogQ]=IndKGStar_TangentPoint(mu, beliefvar, noisevar,nleft)
    M = length(mu);
    assert(length(beliefvar)==M);
    assert(length(noisevar)==M || length(noisevar)==1);
    % We would only measure a perfectly known alternative if all other
    % alternatives were perfectly known, so we can treat "everything perfectly
    % known" as a special-case and otherwise safely ignore those perfectly
    % known alternatives.
    imperfect = find(beliefvar>0);
    if (length(imperfect) == 0)
	disp(sprintf('IndKGStar has a tie, because every alternative is perfectly known.'));
	xkg = ChooseRandomElement(1:M);
	maxLogQ = -Inf;
	return;
    end
    delta = Delta(mu);
    delta = delta(imperfect);

    % Make it work with both scalar noise variances and vectors.
    if (length(noisevar)==1)
	nstar = TangentPoint(delta, beliefvar(imperfect), noisevar, 1, nleft);
	v = noisevar ./ nstar;
	sigmatilde = beliefvar(imperfect) ./ sqrt(beliefvar(imperfect) + v);
    else
	nstar = TangentPoint(delta, beliefvar(imperfect), noisevar(imperfect), 1, nleft);
	v = noisevar(imperfect) ./ nstar;
    	sigmatilde = beliefvar(imperfect) ./ sqrt(beliefvar(imperfect) + v);
    end
    z = -abs(delta) ./ sigmatilde;
    logQ = log(sigmatilde') + LogEI(z) - log(nstar');
    % Now nstar, sigmatilde, z, delta, and logQ are indexed over only the imperfect alternatives.

    % Debugging
    %{
    disp(sprintf('mu=%s', mat2str(mu,3)));
    disp(sprintf('nstar=%s', mat2str(nstar,3)));
    disp(sprintf('sigmatilde=%s', mat2str(sigmatilde,3)));
    disp(sprintf('logQ=%s', mat2str(logQ,3)));
    %}

    xkg = imperfect(argmax(logQ)); % Remember to map back to indexing over all the alternatives.
    maxLogQ = max(logQ);
    % Check for ties, if you wish.
    %if (length(find(logQ==maxLogQ))>1)
    if (length(find(logQ==maxLogQ))>2)
	%disp(sprintf('IndKGStar has a tie.  max(logQ)=%g length(find(logQ==maxLogQ))=%d', maxLogQ, length(find(logQ==maxLogQ))));
	%disp(sprintf('IndKGStar tied alternatives: %s', mat2str(imperfect(find(logQ==maxLogQ)))));
    end
end
