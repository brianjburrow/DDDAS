function value = tMAP(theta, gamma, multi_indices)
    value = 0;
    for idx = 1:length(multi_indices(:,1))
        value = value + gamma(idx) * multivariatePolynomial(theta, ...
            multi_indices(idx,:));
    end
end