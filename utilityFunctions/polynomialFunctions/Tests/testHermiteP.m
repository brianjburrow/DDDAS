function testHermiteP()
    figure(1)
    title("Hermite Polynomials (Probabilists)")
    hold on
    xlim([-4, 7])
    ylim([-10, 20])
    xx = linspace(-4, 7);
    for iOrder = 0:5
        yy = hermiteP(iOrder, xx);
        disp(yy)
        plot(xx, yy)
    end
    legend("n = 0", 'n = 1', 'n = 2', 'n = 3',...
        "n = 4", 'n = 5', 'n = 6', 'n = 7', 'n = 8')
end
        