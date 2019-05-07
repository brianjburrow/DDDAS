function x = benchmark_1_ukfStateFunc(x, time)
    %x = x./2 + 25*x./(1 + x.^2) + 8*cos(1.2*time) + normrnd(0, 10^0.5, 1, length(x));
    x = 1.15*x;
    %x = x./2 + 25*x./(1 + x.^2);
end