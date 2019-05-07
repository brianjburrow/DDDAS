function output = hermiteP(order, x)
    %% Implements Hermite polynomials up to the 8^th order
    if order == 0
        output = ones(size(x));
    elseif order == 1
        output = x;
        
    elseif order == 2
        output = x.^2 - 1;
        
    elseif order == 3
        output = x.^3 - 3*x;
        
    elseif order == 4
        output = x.^4 - 6*x.^2 + 3;
        
    elseif order == 5
        output = x.^5 - 10*x.^3 + 15*x;
        
    elseif order == 6
        output = x.^6 - 15*x.^4 + 45*x.^2 - 15;
        
    elseif order == 7
        output = x.^7 - 21*x.^5 + 105*x.^3 - 105*x;
        
    elseif order == 8
        output = x.^8 - 28*x.^6 + 210*x.^4 - 420*x.^2 + 105;
        
    else
        disp("Polynomial order not coded yet")   
    end
end