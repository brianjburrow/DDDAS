% I was wondering whether ExpMaxNorm was concave in sigma.
% I wrote this function to take random intersections through the space of
% possible sigma and plot.
% It is now probably better implemented using the CheckConcavity function.
function [mu,sigma,deltaSigma] = IsExpMaxConvex(numpoints, M, mu, sigma, deltaSigma)
    if (nargin == 2)
        mu = rand(1,M);
        sigma = rand(1,M);
        deltaSigma = rand(1,M);
    end
    t=[0:.1:1];
    for i=[1:length(t)]
        [y(i),e(i)] = ExpMaxNorm(numpoints,mu,sigma+t(i)*deltaSigma);
    end
    errorbar(t,y,e);
end
