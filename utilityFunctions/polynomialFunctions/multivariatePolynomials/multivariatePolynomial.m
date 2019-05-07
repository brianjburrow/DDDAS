function evaluation = multivariatePolynomial(theta, multiIndex)
    evaluation = 1;
    for idx = 1:length(multiIndex(1,:))
        evaluation = evaluation * hermiteP(multiIndex(1,idx), theta(1,idx));
    end
end