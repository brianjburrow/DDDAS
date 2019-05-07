function legendreP(order, X)
    %% legendreP Implements Legendre Polynomials up to 6th order.
    if order == 0
        output = ones(size(X));
    elseif order == 1
        output = X;
    elseif order == 2
        output = 0.5 * (3*X.^2 - 1);
    elseif order == 3
        output = 0.5 * (5*X.^3 - 3*X);
    elseif order == 4
        output = 0.125 * (35 * X.^4 - 30.*X.^2 + 3);
    elseif order == 5
        output = 0.125 * (63 * X.^5 - 70 * X.^3 + 15 * X);
    else
        error("Legendre Polynomials of 6th order or higher not implemented")
    end
end