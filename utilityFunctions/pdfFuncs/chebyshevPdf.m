function p = chebyshevPdf(X)
    p = 1./(pi.*sqrt(1 - X.^2));
end