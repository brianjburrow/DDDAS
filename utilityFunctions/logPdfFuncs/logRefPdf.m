function p = logRefPdf(r)
    % Compute the log pdf of a multivariate standard normal
    d = length(r);
    
    p = -0.5*d*log(2*pi) - 0.5*sum(r.^2);

end