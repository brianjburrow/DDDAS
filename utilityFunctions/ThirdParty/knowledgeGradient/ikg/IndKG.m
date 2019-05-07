% Implements the independent normal knowledge gradient policy.
% mu should be an M-vector,
% beliefvar should be an M vector, containing belief variances.
% noisevar should either be a scalar homogeneous noise variance, or an M-vector of the noise variances.
function [xkg,maxLogQ]=IndKG(mu, beliefvar, noisevar)
    M = length(mu);
    mu = ConvertToColumnVector(mu);
    beliefvar = ConvertToColumnVector(beliefvar);
    noisevar = ConvertToColumnVector(noisevar);
    assert(length(beliefvar)==M);
    assert(length(noisevar)==M || length(noisevar)==1);
    % We would only measure a perfectly known alternative if all other
    % alternatives were perfectly known, so we can treat "everything perfectly
    % known" as a special-case and otherwise safely ignore those perfectly
    % known alternatives.
    imperfect = find(beliefvar>0);
    if (length(imperfect) == 0)
	disp(sprintf('IndKG has a tie, because every alternative is perfectly known.'));
	xkg = ChooseRandomElement(1:M);
	maxLogQ = -Inf;
	return;
    end
  
    % Make it work with both scalar noise variances and vectors.
    if (length(noisevar)==1)
    	sigmatilde = beliefvar(imperfect) ./ sqrt(beliefvar(imperfect) + noisevar);
    else
    	sigmatilde = beliefvar(imperfect) ./ sqrt(beliefvar(imperfect) + noisevar(imperfect));
    end
    % The two lines below replace a previous line, which read
    %   delta = Delta(mu(imperfect)); 
    % This had a bug in it: if the best alternative was perfect, then that
    % alternative was left out of the "imperfect" array, and thus
    % Delta(mu(imperfect)) gave the wrong answer.
    delta = Delta(mu);
    delta = delta(imperfect);
    z = -abs(delta) ./ sigmatilde;
    logQ = log(sigmatilde') + LogEI(z);
    % Now sigmatilde, z, delta, and logQ are indexed over only the imperfect alternatives.

    %{
    % Debugging
    disp(sprintf('mu=%s', mat2str(mu,3)));
    disp(sprintf('sigmatilde=%s', mat2str(sigmatilde,3)));
    disp(sprintf('logQ=%s', mat2str(logQ,3)));
    %}

    xkg = imperfect(argmax(logQ)); % Remember to map back to indexing over all the alternatives.
    maxLogQ = max(logQ);
    % Check for ties, if you wish.
    %if (length(find(logQ==maxLogQ))>1)
    if (length(find(logQ==maxLogQ))>2)
	%disp(sprintf('IndKG has a tie.  max(logQ)=%g length(find(logQ==maxLogQ))=%d', maxLogQ, length(find(logQ==maxLogQ))));
	%disp(sprintf('IndKG tied alternatives: %s', mat2str(imperfect(find(logQ==maxLogQ)))));
    end
end
