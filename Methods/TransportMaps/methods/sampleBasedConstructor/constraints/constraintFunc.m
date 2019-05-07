function [value, eqCon] = constraintFunc(gamma, samples, multi_indices, index, dmin)
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
    
    pmvp = @(xx, yy) partialMVP(xx, yy, index);
    
    for idx = 1:K
        for dmx = 1:M
            G(idx, dmx) = pmvp(samples(idx,1:dim),multi_indices(dmx,1:dim));
        end
    end
    
    value = dmin*ones([K,1]) - G*gamma';

end