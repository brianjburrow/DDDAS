function value = tMAP_partial_vectorized(thetas, gamma, multi_indices, index)
    % Partial w.r.t a sample dimension, index
    [~, dim] = size(thetas);
    value = zeros(length(thetas(:,1)), 1);
    pmvp = @(xx) partialMVP_vectorized(thetas, xx, index);
    for idx = 1:length(multi_indices(:,1))
        value = value + gamma(idx) * pmvp(multi_indices(idx,1:dim));
    end
end