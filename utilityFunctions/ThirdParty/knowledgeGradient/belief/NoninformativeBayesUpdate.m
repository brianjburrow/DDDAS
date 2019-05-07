% [mu1,A1,B1] = NoninformativeBayesUpdate(mu0,A0,B0,x0,lambdax,yhat1)
% Created on 8/15/2007 in support of the correlated normal
% knowledge-gradient paper.
% This function does a Bayesian update from a noninformative prior
% \theta \sim \Ncal(\mu0,\Sigma0), where \Sigma0 = A0 + k*B0, k->infty.
% After the measurement of point x0, which gave value yhat1 with
% measurement variance lambdax,
% \theta \sim \Ncal(\mu1,\Sigma1), where \Sigma1 = A1 + k*B1, k->infty.
%
% mu0 should be a column vector.
% A0 and B0 should be symmetric square matrices.
% lambdax and yhat1 should be scalars.
function [mu1,A1,B1] = NoninformativeBayesUpdate(mu0,A0,B0,x0,lambdax,yhat1)
    if (B0(x0,x0) ~= 0)
        mu1 = mu0 + (yhat1 - mu0(x0))*B0(:,x0)/B0(x0,x0);
        A1 = A0 + ((A0(x0,x0)+lambdax)*B0 - A0(:,x0)*B0(x0,:) - B0(:,x0)*A0(x0,:))/B0(x0,x0);
        B1 = B0 - B0(:,x0)*B0(x0,:)/B0(x0,x0);
    else
        mu1 = mu0 + (yhat1 - mu0(x0))*A0(:,x0)/(A0(x0,x0)+lambdax);
        A1 = A0 - A0(:,x0)*A0(x0,:)/(A0(x0,x0) + lambdax);
        B1 = B0;
    end