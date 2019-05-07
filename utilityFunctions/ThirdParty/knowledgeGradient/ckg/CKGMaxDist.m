% Plots max \mu^{n+1}_i and argmax \mu^{n+1} as functions of z, where we
% measure at a fixed x.  Useful for understanding why the KG factor is what it
% is.
function [max_val, max_loc, logq] = CKGMaxDist(mu, Sigma, noisevar, x)
    M = length(mu);
    %disp(sprintf('mu=%s', mat2str(mu,3)));
    %disp(sprintf('diag(Sigma)=%s', mat2str(diag(Sigma),3)));

    z=[-10:.1:10];
    denom = sqrt(Sigma(x,x)+noisevar(x));
    if (denom == 0)
	sigmatilde = zeros(1,M);
    else
	sigmatilde = Sigma(:,x) / denom;
    end
    max_val = zeros(size(z));
    max_loc = zeros(size(z));
    for i=1:length(z)
	[max_val(i), max_loc(i)] = max(mu+z(i)*sigmatilde);
    end
    logq  = LogEmaxAffine(mu,sigmatilde);

    subplot(3,1,1)
    stderr=sqrt(diag(Sigma));
    plot(1:M,mu,'-r',1:M,mu+2*stderr,'--r',1:M,mu-2*stderr,'--r',x,mu(x),'xr')
    xlabel('x');
    ylabel('belief +/- 2*stderr');

    subplot(3,1,2)
    plot(z,max_val)
    xlabel('z');
    ylabel('max(mu_i + z*sigmatilde_i)');

    subplot(3,1,3)
    plot(z,max_loc)
    xlabel('z');
    ylabel('argmax(mu_i + z*sigmatilde_i)');
end
