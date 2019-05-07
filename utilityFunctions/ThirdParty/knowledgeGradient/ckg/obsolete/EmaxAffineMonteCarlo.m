% Evaluates the same thing as AffineEmax, but using monte carlo.  This was
% written as a way to test AffineEmax.
%
function [estimate, error] = EmaxAffineMonteCarlo(a,b,nsamples)
    z = randn(1,nsamples);
    result = zeros(1,nsamples); % preallocate
    for n=[1:nsamples]
        result(n) = max(a + b*z(n));
    end
    estimate = mean(result);
    error = std(result)/sqrt(nsamples-1);
end
