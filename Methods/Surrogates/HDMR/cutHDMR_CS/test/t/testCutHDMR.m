order = 2;
cutLine = [0.5, 0.5];
nSamples = 10;
fm = @(xx) fullModel(xx);

hdmrModel = cutHDMR(fm, order, cutLine, nSamples);

function point = fullModel(X)
    point = sum(X, 2);
end