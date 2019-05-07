function data = cappe_ukfDataGen_noisy(x, time)
    data = 0.05 * x.^2 + normrnd(0, 1, size(x));
end