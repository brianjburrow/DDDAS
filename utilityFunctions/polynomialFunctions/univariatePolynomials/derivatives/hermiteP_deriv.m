function output = hermiteP_deriv(order, x)
    if order == 0
        output = 0;
    elseif order == 1
        output = 1;
    elseif order == 2
        output = 2*x;
    elseif order == 3
        output = 3*x.^2 - 3;
    elseif order == 4
        output = 4*x.^3 - 12*x;
    elseif order == 5
        output = 5*x.^4 - 30*x.^2 + 15;
    elseif order == 6
        output = 6*x.^5 - 15*4*x.^3 + 45*2 * x;
    elseif order == 7
        output = 7*x.^6 - 21*5*x.^4 + 105*3*x.^2 - 105;
    elseif order == 8
        output = 8*x.^7 - 28*6*x.^5 + 210*4 * x.^3 - 420*2*x;
    else
        disp("Polynomial order not coded yet")   
    end
end