clear all
close all
clc
plotUB = 100;
numMeas = 100;

p = sobolset(3);

MU = p(1:numMeas,:);

MU(:,1) = forwardModel(120*MU(:,2) - 60);
MU(:,2) = 120*MU(:,2) - 60;
MU(:,3) = 60^2*MU(:,3) + 100;

burnIn = 5000;
chainLength = 10000;
samples = zeros(chainLength, numMeas);
counter = 1;

for dmx = 1:numMeas
    fprintf("Iteration Number %d \n", dmx)
    samples(1, dmx) = 0;
    C = 60^2;
    measurement = MU(dmx, :);
    targetPdf = @(xx) likelihood(measurement, xx);
    for idx = 2:chainLength
        proposalState = normrnd(samples(idx-1, dmx), C);
        propPdf = @(xx, yy, zz) logmvnpdf(xx, samples(idx-1, dmx), zz);
        [samples(idx, dmx), accept, C] = metropolisHastingsIter(...
            samples(idx-1, dmx), ...
            proposalState,...
            propPdf, ...
            targetPdf, ...
            idx, ...
            C, ...
            samples(1:idx, dmx));
    end
end

samples = samples(burnIn:end, :);

for idx = 1:length(samples(1,:))
    if idx >1
        cla
    end
    hist(samples(:,idx))
    xlim([-65, 65])
    drawnow
    pause(1)
end
    

function p = likelihood(data, state)
    mu = forwardModel(state);
    p  = logmvnpdf(data(1,1), mu, 1) + logmvnpdf(data(1,1), data(1,2), data(1,3)); %+ logmvnpdf(state, 2, 1);
end

function meas = forwardModel(input)
    %meas = 1*input + sin(input);
    meas = 0.05*input.^2;
end
