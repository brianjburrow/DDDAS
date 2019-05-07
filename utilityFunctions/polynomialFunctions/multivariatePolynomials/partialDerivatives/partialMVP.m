function evaluation = partialMVP(theta, multiIndex, index)
    n = multiIndex(1,index);                                                  % Index is the variable we take the derivative w.r.t
    
    evaluation = hermiteP_deriv(n, theta(index));

    for idx = 1:length(multiIndex(1,:))
        if idx ~= index
            evaluation = evaluation * hermiteP(multiIndex(1,idx), theta(1,idx));  
        end
    end
end