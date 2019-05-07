% Implements the knowledge gradient policy for the correlated normal problem.
% Uses the covariance matrix sigma rather than its square root, or the
% precision matrix beta.
% mu should be an M-vector,
% beta should be an MxM matrix.
% noisevar should be an M-vector of the variances.
function [xkg,maxLogQ,logQ]=CorrelatedNormalKGWithSigma(mu, Sigma, noisevar)
    M = length(mu);
    %disp(sprintf('mu=%s', mat2str(mu,3)));
    %disp(sprintf('diag(Sigma)=%s', mat2str(diag(Sigma),3)));
    logQ = zeros(1, M);
    for x = 1:M
        denom = sqrt(Sigma(x,x)+noisevar(x));
        if (denom == 0)
            logQ(x) = -Inf;
        else
            sigmatilde = Sigma(:,x) / denom;
            % Same as: 
            % sigmatilde = Sigma*e(x,M) / sqrt(Sigma(x,x)+noisevar(x));
                %Q(x) = EmaxAffine(mu,sigmatilde);
            logQ(x) = LogEmaxAffine(mu,sigmatilde);
            %disp(sprintf('x=%d logQ=%f sigmatilde=%s', x, logQ(x), mat2str(sigmatilde,3)));
        end
    end
    xkg = Argmax(logQ);
    maxLogQ = max(logQ);
    % Check for ties, if you wish.
    if (length(find(logQ==maxLogQ))>1)
        %warning('Correlated KG has a tie for alternative with the largest Q-factor.  Choosing uniformly at random among them.');
	fprintf('CKG has a tie.  max(logQ)=%g length(find(logQ==maxLogQ))=%d\n', maxLogQ, length(find(logQ==maxLogQ)));
	fprintf('CKG tied alternatives: %s\n', mat2str(find(logQ==maxLogQ)));
	%disp(Sigma);
    end
end
