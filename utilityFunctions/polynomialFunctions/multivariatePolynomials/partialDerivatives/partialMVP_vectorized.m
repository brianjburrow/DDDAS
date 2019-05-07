function evaluation = partialMVP_vectorized(thetas, multiIndex, index)
    n = multiIndex(1,index);                                                % Index is the variable we take the derivative w.r.t
    
    evaluation = hermiteP_deriv(n, thetas(:,index));

    
    for idx = 1:length(multiIndex(1,:))
        if idx ~= index
            evaluation = evaluation .* hermiteP(...
                multiIndex(1,idx),...
                thetas(:,idx)...
                );  
        end
    end
end