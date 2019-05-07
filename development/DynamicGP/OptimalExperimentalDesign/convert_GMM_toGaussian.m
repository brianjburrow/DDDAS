function [mean, std_dev] = convert_GMM_toGaussian(XX, means, sigmas)
    % XX is dummy variable, don't want to remove or it may break old
    % codes
    num_dists = length(means);
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