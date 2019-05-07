function testBesselP()
    figure(1)
    title("Bessel Polynomials")
    hold on
    xlim([-2, 2])
    ylim([-6, 6])
    xx = linspace(-2, 2);
    for iOrder = 0:5
        yy = besselP(iOrder, xx);
        disp(yy)
        plot(xx, yy)
    end
    legend("n = 0", 'n = 1', 'n = 2', 'n = 3',...
        "n = 4", 'n = 5')
end
        