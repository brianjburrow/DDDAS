% Calculates the standard deviation of the posterior mean that results from
% measuring an alternative with a given prior variance and given noise
% variance.  Note that this returns 0 if noiseVariance is infinite.
function sigmatilde = SigmaTilde(priorVariance, noiseVariance)
	sigmatilde = sqrt(priorVariance ./ (1 + noiseVariance./priorVariance));
end
