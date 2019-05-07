function x = cappe_ukfStateFunc(x, time)
    %x = x./2 + 25*x./(1 + x.^2) + 8*cos(1.2*time) + normrnd(0, 10^0.5, 1, length(x));
    x = x./2 + 25*x./(1 + x.^2) + 8*cos(1.2*time);
    %x = x./2 + 25*x./(1 + x.^2);
end