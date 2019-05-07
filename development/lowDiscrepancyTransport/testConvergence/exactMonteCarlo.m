function theta = random_sampler(numSamples)
    %theta = zeros([2,numSamples]);
    r_1 = normrnd(0, 1, [1, numSamples]);
    r_2 = normrnd(0, 1, [1, numSamples]);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end

function output = forwardModel()

end