function multi_indices = genNoCrossMI(order, dimensions)
    % order is the polynomial order that we want to use
    % dimensions is the dimension of the problem
    % limit limits the number of columns that we fill
    % if limit > dimensions, then we fill all columns

    for idx = 1:dimensions
        initVec(idx) = order+1;
    end
    options = fullfact(initVec);
    
    options = options - ones(size(options));
    
    count = 1;
    for idx = 1:length(options(:,1))
        bincheck = checkSetInclusion(options(idx,:), order, dimensions);
        if bincheck
            multi_indices(count,:) = options(idx,:);
            count = count + 1;
        end
    end
end

function binaryCheck = checkSetInclusion(VECTOR, order, dim)
    magSum = norm(VECTOR, 1);
    if (magSum <= order) && (magSum >= 0)
        binaryCheck = true;
        for idx = 1:length(VECTOR(1,:))
            if  (VECTOR(1,idx) ~= 0) && (binaryCheck) && (idx ~= dim)
                binaryCheck = false;
            end
        end
    else
        binaryCheck = false;
    end
end