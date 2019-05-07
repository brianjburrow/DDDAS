function output = besselP(order, x)
    if order == 0
        output = ones(size(x));
    elseif order == 1
        output = x + 1;
    elseif order == 2
        output = 3*x.^2 + 3*x + 1;
    elseif order == 3
        output = 15*x.^3 + 15*x.^2 + 6*x + 1;
    elseif order == 4
        output = 105*x.^4 + 105*x.^3 + 45*x.^2 + 10*x + 1;
    elseif order == 5
        output = 945 * x.^5 + 945*x.^4 + 420*x.^3 + 105 * x.^2 + 15*x + 1;
    else
        disp("Higher order Bessel olynomials not coded yet")   
    end
end