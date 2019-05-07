function samples = genPosteriorSamples(measurements, numSamples, burnIn, thining)

slope = 0.5;
additive = 5;
dataSlope = 3;
plotUB = 100;
numMeas = length(measurements);

xx = linspace(0, 1, 1000);
mu = measurements;

chainLength = numSamples + burnIn-1;
samples = zeros(chainLength, numMeas);
counter = 1;
for dmx = 1:numMeas
    samples(1, dmx) = 0;
    C = 0.5;
    measurement = mu(dmx);
    targetPdf = @(xx) log_likelihood(measurement, xx);
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

samples = samples(burnIn:thining:end, :);

function p = log_likelihood(data, state)
    mu = forwardModel(state);
    p  = logmvnpdf(data, mu, 1); %+ logmvnpdf(state, 2, 1);
end

function meas = forwardModel(input)
    %meas = 1*input + sin(input);
    meas = input^2;
end
end