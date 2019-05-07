function [value, eqCon] = constraint_deriv(gamma, samples, multi_indices, index, ~)
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

    value_dir = zeros([dim, 1]);
    vec = zeros(size(gamma));
    for idx = 1:dim
        vec(idx) = 1;
        value_dir(idx) = -G*vec';
        vec(idx) = 0;
    end
end