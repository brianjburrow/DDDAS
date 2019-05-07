clear all 
close all
nSamples = 25;
lowBound = -1.5;
upBound = 1.5;
xx = linspace(-5,10, nSamples);
yy = linspace(0,15, nSamples);
[XX, YY] = meshgrid(xx, yy);
zz = zeros(size(XX));
for iColumn = 1:length(zz(1,:))
    zz(:, iColumn) = BraninFunction(XX(:,iColumn), YY(:, iColumn));
end
surf(XX, YY, zz)
hold on
title("Top down view of the surface plot")
colorbar()
view(2)

figure(2)
xlim([-5, 10])
hold on
ylim([0, 15])
title("Quiver Plot of the Derivative")
gradZZ1 = zeros(size(XX));
gradZZ2 = zeros(size(XX));
for iColumn = 1:length(zz(1,:))
    [gradZZ1(:, iColumn), gradZZ2(:, iColumn)] = BraninGradient(XX(:,iColumn), YY(:, iColumn));
end


x = reshape(XX, nSamples^2, 1);
y = reshape(YY, nSamples^2, 1);
gradzz1 = reshape(gradZZ1, nSamples^2, 1);
gradzz2 = reshape(gradZZ2, nSamples^2, 1);

quiver(x, y, gradzz1, gradzz2)

function output = BraninFunction(x1, x2)
    % https://www.sfu.ca/~ssurjano/branin.html
    a = 1;
    b = 5.1 / (4*pi^2);
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    [nSamples, ~] = size(x1);
    oneVec = ones(nSamples, 1);
    output = a * (x2 - b * x1.^2 + c*x1 - r*oneVec).^2 + ...
        s*(1 - t)*cos(x1) + s*oneVec;
end

function [part1, part2] = BraninGradient(x1, x2)
    a = 1;
    b = 5.1 / (4*pi^2);
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    [nSamples, ~] = size(x1);
    oneVec = ones(nSamples, 1);
    % partial calculation
    repeatedCalc = (x2 - b*x1.^2 + c*x1 - r*oneVec);
    % component 1 
    part1 = (2*a*repeatedCalc.*(-2*b*x1 + c*oneVec) - s*(1 - t)*sin(x1))';
    part2 = (2*a*repeatedCalc)';

end