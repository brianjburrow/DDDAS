clear all
close all
clc

% Samples will be nSamples x nDimensions

% Results stored as trial(iTrial).samples(iTime, iRepetition)
% where iTrial is the current trial, which controls parameter values
% iTime is the current time step in the dynamical system
% and iRepetition is the current repetition handling the stochastic evolution

nTrials = 50; nTime = 30; nRepetitions = 50;

trial(10).samples = zeros(nTime, nRepetitions);

for iTrial = 1:nTrials
    trial(iTrial).samples = zeros(nTime, nRepetitions);
    params = hyperSamples(1);                                               % Generate parameters for this trial   
    trial(iTrial).samples(1,:) = normrnd(...
        0,(1/(1 - params(1,2)^2))^0.5, [1, nRepetitions]...                 % initialize the repetitions for each trial
        );
    for iTime = 2:nTime
        trial(iTrial).samples(iTime, :) = evolve(...
            trial(iTrial).samples(iTime - 1, :), params...                  % Perform the state evolution for each repetition simultaneously
            );
    end
end

samples = unpackWhole(trial);
[a,b] = size(samples);
samples = reshape(samples, [1, a*b]);
maxVar = var(samples);
maxStd = std(samples);

minStd = 1;
raisedSamples = normrnd(samples(1), maxStd, [1,100]);
disp('Computing raised samples')
inSamp = 100;
raisedSamples = zeros(1, length(samples)*inSamp);
for iSamp = 1:length(samples)
    raisedSamples(((iSamp-1)*inSamp + 1):(iSamp*inSamp)) = normrnd(samples(iSamp), maxStd, [1,inSamp]);
end
disp('Completed Raised samples')

bounds = zeros(1,2);
bounds(1,1) = min(raisedSamples);
bounds(1,2) = max(raisedSamples);

xx = linspace(min(samples), max(samples), 1000);
figure(1)
histogram(samples, 'Normalization', 'pdf')
hold on
plot(xx, normpdf(xx, 0, maxStd), 'red', 'LineWidth', 3.0)
hold on
plot(xx, mixturePdf(xx, samples, maxStd*ones(size(samples))), 'blue', 'LineWidth', 3.0)
xlabel('State')
ylabel('Probability Density')

% optimize transport map

ppSet = mixturePdf(raisedSamples(:,1:1000:end), samples, maxStd*ones(size(samples)));
disp(min(ppSet))
threshold = 40;

% nBase2 = determineSampleSizeWeighted(samples', maxStd*ones(size(samples')), threshold, bounds, min(ppSet));
% 
% nBase3 = determineSampleSizeWeighted(samples', ones(size(samples')), threshold, bounds, min(ppSet));

nBase = determineSampleSize(samples', maxStd*ones(size(samples')), threshold, bounds);

%% Gaussian Map
propMu = mean(raisedSamples);
propStd = std(raisedSamples);

nGauss = determineGaussianSampleSize(propMu, propStd, samples', maxStd*ones(size(samples')), threshold, nBase);

for dmx = 1:numDim
    fprintf("DMX: %d \n", dmx)
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder,  dmx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder,  dmx);
    else
        multi_indices = genNoCrossMI(polyOrder,  dmx);
    end
    
    MD = length(multi_indices(:,1));
    value = min(find(multi_indices(:,dmx) > 0));
    init = zeros(1, MD);
    init(value) = 1;
    %init(2) = 1;
    %init = [ones([1, value]), zeros([1, MD-value])];
    if dmx > 1
        numGamPrev = length(gammas(dmx-1).val);
        init(1:numGamPrev) = gammas(dmx-1).val;
    end
    gammas(dmx).val = optimizeMap_ESS(samples(:,1:dmx),...
                                     [],...
                                     multi_indices, init);
end
    refSamples = zeros(size(samples));
    tic
for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder,  idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder,  idx);
    else
        multi_indices = genNoCrossMI(polyOrder,  idx);
    end
    refSamples(:,idx) = tMAP_vectorized(samples(:,1:idx), gammas(idx).val, multi_indices); 
end


function n = determineSampleSizeWeighted(muSet, stdSet, threshold, bounds, pThresh)
    nDim = length(muSet(1,:));
    nSet = length(muSet(:,1));
    ss = sobolset(nDim);
    n = 2;
    cont = 1;
    maxim = bounds(1,2);
    minim = bounds(1,1);
    gap = maxim - minim;  
    samples = ss(1, :)*gap + minim;
    while cont
        cont2 = 1;
        while cont2
            newSamp = ss(n, :)*gap + minim;
            pp = mixturePdf(newSamp, muSet, stdSet);
            if pp > pThresh
                samples = [samples; newSamp];
                cont2 = 0;
            else
                n = n + 1;
            end
        end
        fprintf('max: %f, min: %f \n', max(samples), min(samples))
        arrNEff = zeros(1,nSet);
        for iSet = 1:nSet
            arrNEff(iSet) = computeUniformToNormalESS(samples,...
                muSet(iSet, :),...
                stdSet(iSet, :));
        end
        worstCase = min(arrNEff);
        fprintf('Worst Case EFF = %f \n', worstCase);
        bestCase = max(arrNEff);
        fprintf('Best Case EFF = %f \n', bestCase);
        if worstCase >= threshold
            cont = 0;                                                       % break while loop
        else
            n = n + 1;
        end
        fprintf('Continuing, new n = %d \n', length(samples(:,1)))
        fprintf('True n = %d \n', n)
    end
end

function n = determineSampleSize(muSet, stdSet, threshold, bounds)
    nDim = length(muSet(1,:));
    nSet = length(muSet(:,1));
    ss = sobolset(nDim);
    n = 10*nDim;
    cont = 1;
    maxim = bounds(1,2);
    minim = bounds(1,1);
    gap = maxim - minim;
    
    while cont
        samples = ss(1:n, :)*(gap) + minim;
        fprintf('max: %f, min: %f \n', max(samples), min(samples))
        arrNEff = zeros(1,nSet);
        for iSet = 1:nSet
            arrNEff(iSet) = computeUniformToNormalESS(samples,...
                muSet(iSet, :),...
                stdSet(iSet, :));
        end
        worstCase = min(arrNEff);
        fprintf('Worst Case EFF = %f \n', worstCase);
        bestCase = max(arrNEff);
        fprintf('Best Case EFF = %f \n', bestCase);
        if worstCase >= threshold
            cont = 0;                                                       % break while loop
        else
            n = n + 10;
        end
        fprintf('Continuing, new n = %d \n', n)
    end
end

function nEff = computeUniformToNormalESS(samples, mu, std)
    nSamp = length(samples(:,1));
    weights = normpdf(samples, mu, std)';
    weights = weights./nSamp;
    sWeights = sum(weights)^2;
    sWeightSq = sum(weights.^2);
    nEff = sWeights/sWeightSq;
end

function n = determineGaussianSampleSize(propMu, propStd, muSet, stdSet, threshold, nMax)
    nDim = length(muSet(1,:));
    nSet = length(muSet(:,1));
    ss = sobolset(nDim);
    n = nMax;
    cont = 1;
    
    while cont
        samples = uniformTOgaussian(ss(1:n));
        samples = samples(1:n, :)*(propStd) + ones(n, nDim)*propMu;
        fprintf('max: %f, min: %f \n', max(samples), min(samples))

        arrNEff = zeros(1,nSet);
        for iSet = 1:nSet
            arrNEff(iSet) = computeUniformToNormalESS(samples,...
                muSet(iSet, :),...
                stdSet(iSet, :));
        end
        worstCase = min(arrNEff);
        fprintf('Gauss Worst Case EFF = %f \n', worstCase);
        bestCase = max(arrNEff);
        fprintf('Gauss Best Case EFF = %f \n', bestCase);
        if worstCase > threshold
            cont = 0;                                                       % break while loop
        else
            n = n + 10;
            propStd = propStd + 1;
        end
        fprintf('Continuing, new n = %d \n', n)
    end
end

function nEff = computeGaussianToNormalESS(samples, propMu, propStd, mu, std)
    nSamp = length(samples(:,1));
    weights = normpdf(samples, mu, std)';
    weights = weights./(normpdf(samples, propMu, propStd)');
    sWeights = sum(weights)^2;
    sWeightSq = sum(weights.^2);
    nEff = sWeights/sWeightSq;
end

function p = mixturePdf(xx, muSet, stdSet)
    p = normpdf(xx, muSet', stdSet');
%     for iDist = 1:length(muSet)
%         p = p + normpdf(xx, muSet(iDist), stdSet(iDist));
%     end
    p = sum(p)./length(muSet);
end

function samps = unpackWhole(trialObj)
    nTrials = length(trialObj);
    samps = trialObj(1).samples;
    for iTrial = 2:nTrials
        samps = [samps, trialObj(iTrial).samples];
    end
end

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end

function inits = hyperSamples(nSamp)
    sigM = 0.01;
    mu = normrnd(0, sigM, [nSamp,1]);
    phi = betarnd(10, 1, [nSamp,1]).*2 - ones([nSamp, 1]);
    inits = [mu, phi];
end

function state = evolve(state, params)
    mu = params(1,1);
    phi = params(1,2);
    state = mu + phi*(state - mu) + normrnd(0,1, size(state));
end

function measurement = forwardModel(state)
    %measurement = normrnd(0,1, size(state)).*exp(state/2.0);
    measurement = exp(state/2.0);
end
