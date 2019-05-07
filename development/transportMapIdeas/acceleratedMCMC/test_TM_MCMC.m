clear all
clc
close all

rng(0)
polyOrder = 3;

measurements = 2;
numInitSamples = 20000;
numTotalSamples = 30000;
burnIn = 5000;
thining = 5; 

numMeas = length(measurements);

xx = linspace(0, 1, 1000);
mu = measurements;

adaptSpacing = 100;

chainLength = numInitSamples + burnIn-1;
samples = zeros(chainLength, numMeas);
counter = 1;
for dmx = 1:numMeas
    samples(1, dmx) = 0;
    C = 0.5;
    measurement = mu(dmx);
    targetPdf = @(xx) log_likelihood(measurement, xx) + log_prior(xx);
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

initSample = samples(burnIn:thining:end, :);

log_like = @(xx) log_likelihood(measurement, xx);
log_pri = @(xx) log_prior(xx);

disp('outside')
tm_sampler = transportMapMCMC(...
    log_like,...
    log_pri,...
    initSample,...
    polyOrder,...
    adaptSpacing,...
    numTotalSamples...
    );

function p = log_likelihood(data, state)
    mu = forwardModel(state);
    p  = logmvnpdf(data, mu, 1); %+ logmvnpdf(state, 2, 1);
end

function log_p = log_prior(state)
    log_p = logmvnpdf(state, 0, 10);
end

function meas = forwardModel(input)
    %meas = 1*input + sin(input);
    meas = input^2;
end