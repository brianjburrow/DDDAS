clear all
close all

numRepeats = 100;
sampleSizes = 10:10:300;

truth = mean(forwardModel(random_sampler(10000000)));

mean1 = zeros(numRepeats, 2);
for iRepeat = 1:numRepeats
    for iSampSize = 1:length(sampleSizes)
        samples1 = random_sampler(sampleSizes(iSampSize));
        outputs = forwardModel(samples1);
        mean1(iRepeat,iSampSize) = abs(mean(outputs) - truth);
    end
end

std1 = std(mean1);
mean1 = mean(mean1);


mean3 = zeros(numRepeats, 2);
for iRepeat = 1:numRepeats
    for iSampSize = 1:length(sampleSizes)
        samples2 = random_sampler_3(sampleSizes(iSampSize));
        outputs = forwardModel(samples2);
        mean3(iRepeat,iSampSize) = abs(mean(outputs) - truth);
    end
end
std3 = std(mean3);
mean3 = mean(mean3);

figure(1)
%subplot(3,1,1)
errorbar(sampleSizes, mean1, 3*std1, 'red', 'LineWidth', 2.0)
hold on
%subplot(3,1,3)
errorbar(sampleSizes, mean3, 3*std3, 'black', 'LineWidth', 2.0)
hold on
ylabel("Absolute Error")
xlabel("Sample Size")
legend("Target Random", "Sobol Transported")
hold on

figure(2)
subplot(2,1,1)
scatter(samples1(:,1), samples1(:,2),10, 'filled')
hold on
xlim([1, 9])
ylim([1, 9])
xlabel("X_1")
ylabel("X_2")
subplot(2,1,2)
scatter(samples2(:,1), samples2(:,2),10, 'filled')
hold on
xlim([1, 9])
ylim([1, 9])
xlabel("X_1")
ylabel("X_2")


%% Test LDS through Approximate Map
% Transport Map Parameters
rng(0)
polyOrder = 3;
multi_index_type = "TO";                                                    % "TO" Total order, "NM" No Mixed, "NC" No cross
numSamples = 5000;

samples = random_sampler(numSamples);

[~, numDim] = size(samples);

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
    gammas(dmx).val = optimizeMap(samples(:,1:dmx),...
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
fprintf("Computing Reference Samples: %d \n",toc)

tic
for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder,  idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder,  idx);
    else
        multi_indices = genNoCrossMI(polyOrder,  idx);
    end
    betas(idx).val = invertMap(samples(:,1:idx), refSamples(:,1:idx), gammas(idx).val, multi_indices);
end
fprintf("Constructing Inverse Map: %d \n",toc)

tic
mean4 = zeros(numRepeats, 2);
for iRepeat = 1:numRepeats
    for iSampSize = 1:length(sampleSizes)
        ss = sobolset(2);
        ss = scramble(ss, 'MatousekAffineOwen');
        uniSamples = ss(1:sampleSizes(iSampSize), :);                       % Generate Quasi-random uniform samples
        refSamples = uniformTOgaussian(uniSamples);                         % Map them to Quasi-random Gaussian samples
        SAMPLES = zeros(sampleSizes(iSampSize), 2);
        for idx = 1:numDim
            if multi_index_type == "TO"
                multi_indices = genTotalOrderMI(polyOrder, idx);
            elseif multi_index_type == "NM"
                multi_indices = genNoMixedMI(polyOrder, idx);
            else
                multi_indices = genNoCrossMI(polyOrder, idx);
            end
            SAMPLES(:,idx) = tMAP_vectorized(...
                refSamples(:,1:idx),...
                betas(idx).val,...
                multi_indices...
                ); 
        end
        outputs = forwardModel(SAMPLES);
        mean4(iRepeat, iSampSize) = abs(mean(outputs) - truth);
    end
end
std4 = std(mean4);
mean4 = mean(mean4);


figure(1)
%subplot(3,1,1)
errorbar(sampleSizes, mean4, 3*std4, 'blue', 'LineWidth', 2.0)
legend("Target Random", "Exact Sobol Transport", "Approximate Sobol Transport")

%% Begin Function Definitions
function output = forwardModel(x)
    output = x(:,1).^2 + cos(x(:,2));
end

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end

function theta = random_sampler(numSamples)
    %theta = zeros([2,numSamples]);
    r_1 = normrnd(0, 1, [1, numSamples]);
    r_2 = normrnd(0, 1, [1, numSamples]);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end


function theta = random_sampler_3(numSamples)
    %theta = zeros([2,numSamples]);
    ss = sobolset(2);
    ss = scramble(ss, 'MatousekAffineOwen');
    uniSamples = ss(1:numSamples, :);
    refSamples = uniformTOgaussian(uniSamples)';
    r_1 = refSamples(1,:);
    r_2 = refSamples(2,:);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end