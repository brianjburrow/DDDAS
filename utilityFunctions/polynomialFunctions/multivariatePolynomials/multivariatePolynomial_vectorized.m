function evaluation = multivariatePolynomial_vectorized(theta, multiIndex)
    evaluation = ones(length(theta(:,1)), 1);
    for idx = 1:length(multiIndex(1,:))
        evaluation = evaluation .* hermiteP(multiIndex(1,idx), theta(:,idx));
    end
end