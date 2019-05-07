function value = tMAP_vectorized(thetas, gamma, multi_indices)
    [~, dim] = size(thetas);
    value = zeros(length(thetas(:,1)), 1);
    mvp = @(xx) multivariatePolynomial_vectorized(thetas, xx);
    for idx = 1:length(multi_indices(:,1))
        value = value + gamma(idx) * mvp(multi_indices(idx,1:dim));
    end
end