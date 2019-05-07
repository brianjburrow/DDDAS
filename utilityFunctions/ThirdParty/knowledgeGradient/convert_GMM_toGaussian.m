function [mean, std_dev] = convert_GMM_toGaussian(XX, means, sigmas)

    num_dists = length(means);
    value = 0;
    for idx = 1:num_dists
        value = value + normpdf(XX, means(idx), sigmas(idx));
    end
    value = value / num_dists;
    mean = 0;
    variance = 0;
    for idx = 1:num_dists
        mean = mean + means(idx);
    end
    mean = mean/num_dists;
    for idx = 1:num_dists
        variance = variance + (sigmas(idx) + (means(idx) - mean)^2);
    std_dev = (variance / num_dists)^0.5;
    end
end