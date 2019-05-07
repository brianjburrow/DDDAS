function gammas = optimizeMap(samps, prevTerms, multi_indices, init)
    [~, dim] = size(samps);
    
    options = optimoptions(...
                'fmincon', 'MaxIterations',  1000, 'algorithm',...
                'sqp', 'StepTolerance', 1000*eps,...
                'MaxFunctionEvaluations', 10000,...
                'SpecifyObjectiveGradient', true,...
                'SpecifyConstraintGradient',true);
            
    objF  = @(xx) objectiveFunc_vectorized([prevTerms, xx],...
                    samps, multi_indices, dim);
    
    const = @(xx) constraintFunc_vectorized([prevTerms, xx],...
                    samps, multi_indices, dim, 10^-8);                       % Code derivative information into this problem
    
    gammas = fmincon(...
                objF, ...
                init,...                                                     % change "1" to predicted terms eventually
                [],[],[],[],[],[],...
                const, ...
                options...
                );
end