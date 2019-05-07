% [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises)
% The inputs x and y are vectors.
% The input noises can either be a scalar or a vector of the same length as x.
% Corrects newmu(x) to y and newSigma(x) to 0 on perfect measurements.
% SigmaSqrt should be lower triangular.
function [newmu, newSigma] = BayesUpdateBlock(mu, Sigma, x, y, noises)
	% cholinc may produce warning messages if Sigma has 0s on the diagonal,
	% as it will if we are doing perfect measurements.  To avoid these
	% warning messages, do not fix diagonal covariance entries of perfect
	% measurements to 0 until after all the measurements are done.
	% cholinc gives upper triangular, and the transpose switches to lower.
	SigmaSqrt = ichol(sparse(Sigma),1e-15)'; 
	[newmu, newSigmaSqrt]=BayesUpdateSqrtBlock(mu, SigmaSqrt, x, y, noises);
	newSigma = newSigmaSqrt*newSigmaSqrt';

	% Post-procissing fixing of mu and Sigma on perfect measurements.
	perfect = find(noises==0);
	k = length(perfect);
	if (k>0)
		mudiff = max(abs(newmu(x(perfect))-y(perfect)));
		if (mudiff>0)
			disp(sprintf('Fixing mu on perfect measurements, max(mudiff)=%g', mudiff));
			assert(mudiff<1);
			newmu(x(perfect)) = y(perfect);
		end
		sigmadiff = max([max(abs(newSigma(x(perfect),:))),max(abs(newSigma(:,x(perfect))))]);
		if (sigmadiff>0)
			assert(sigmadiff<1);
			disp(sprintf('Fixing Sigma on perfect measurements, max(sigmadiff)=%g', full(sigmadiff)));
            M = length(mu);
            newSigma(x(perfect),:) = zeros(k,M);
			newSigma(:,x(perfect)) = zeros(M,k);
		end
    	end
