function [value, eqCon, DC, DCeq] = constraintFunc_vectorized(gamma, ...
    samples, multi_indices, index, dmin...
    )
    %% Variables
    % Samples:       from target distribution
    % Gamma:         parameters of the transport map
    % Index:         Take derivatives w.r.t
    % Multi_indices: Determines polynomial expansion
    [K, dim] = size(samples);
    M = length(multi_indices(:,1));
    
    eqCon = [];
    %% Construct F,G Matrices
    G    = zeros([K, M]);
    
    pmvp = @(xx, yy) partialMVP_vectorized(xx, yy, index);
    
    for dmx = 1:M
        G(:, dmx) = pmvp(samples(:,1:dim),multi_indices(dmx,1:dim));
    end

    
    value = dmin*ones([K,1]) - G*gamma';
     
    if nargout > 2
        DCeq = [];
        DC = -G';
    end
end