function testLegendreP()
    figure(1)
    title("Legendre Polynomials")
    hold on
    xlim([-1, 1])
    ylim([-1.1, 1.1])
    xx = linspace(-2, 2);
    for iOrder = 0:5
        yy = legendreP(iOrder, xx);
        disp(yy)
        plot(xx, yy)
    end
    legend("n = 0", 'n = 1', 'n = 2', 'n = 3',...
        "n = 4", 'n = 5')
end
       