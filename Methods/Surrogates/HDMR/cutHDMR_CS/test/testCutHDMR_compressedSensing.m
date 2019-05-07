clear all
close all
clc

inputSampler = @(numSamp) [chebyshevSampler(numSamp, -1, 1), chebyshevSampler(numSamp, -1, 1)];
cutLine = [0, 0];
nSamples = 3;
fullModel = @(X) testFunc(X);
hdmr = cutHDMR_compressedSensing(fullModel, cutLine, nSamples, 4,...
    'legendre', inputSampler);
hdmr = hdmr.run();
hdmr.plotSingleComponentFunction(1, 1000)
hdmr.plotSingleComponentFunction(2, 1000)

%% plot test function
% xx = linspace(-1,1);
% yy = linspace(-1,1);
% 
% [XX, YY] = meshgrid(xx, yy);
% 
% ZZ = zeros(size(XX));
% 
% for iX = 1:100
%     for iY = 1:100
%         ZZ(iX, iY) = testFunc([XX(iX, iY), YY(iX, iY)]);
%     end
% end
% 
% surf(XX, YY, ZZ)

function output = testFunc(X)
    output = 1 + X(:,1) + 0.5*(3*X(:,2).^2 - 1) + X(:,1).*X(:,2);
end
