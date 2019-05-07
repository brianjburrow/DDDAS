clear all
close all
% Transport Map Parameters with composite maps
rng(0)
polyOrder = 3;
multi_index_type = "TO";    % "TO" Total order, "NM" No Mixed, "NC" No cross
numSamples = 5000;

samples = random_sampler(numSamples);
%samples = load("Banana_samples.m"); 
%samples = [normrnd(10, 2, [numSamples,1]), normrnd(10,2, [numSamples,1])];

[~, numDim] = size(samples);
figure(1)
scatterhist(samples(:,1), samples(:,2))
hold on
title("True Target")
tic

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
fprintf("Constructing Forward Map: %d \n",toc)

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
figure(2)
scatterhist(refSamples(:,1), refSamples(:,2))
hold on
title("Reference Samples")


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
refSamples = mvnrnd([0,0], eye(2), numSamples);
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
fprintf("Computing Inverted Samples: %d \n",toc)
hold on
figure(3)
scatterhist(SAMPLES(:,1), SAMPLES(:,2))
title("Low Discrepancy Reconstructed Targets")

%% Generating Second MAP
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
        numGamPrev = length(gammas2(dmx-1).val);
        init(1:numGamPrev) = gammas2(dmx-1).val;
    end
    gammas2(dmx).val = optimizeMap(refSamples(:,1:dmx),...
                                     [],...
                                     multi_indices, init);
end
fprintf("Constructing Forward Map: %d \n",toc)

refSamples2 = zeros(size(samples));
tic
for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder,  idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder,  idx);
    else
        multi_indices = genNoCrossMI(polyOrder,  idx);
    end
    refSamples2(:,idx) = tMAP_vectorized(refSamples(:,1:idx), gammas2(idx).val, multi_indices); 
end
fprintf("Computing Reference Samples: %d \n",toc)
figure(2)
scatterhist(refSamples2(:,1), refSamples2(:,2))
hold on
title("Reference Samples")


tic
for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder,  idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder,  idx);
    else
        multi_indices = genNoCrossMI(polyOrder,  idx);
    end
    betas2(idx).val = invertMap(refSamples(:,1:idx), refSamples2(:,1:idx), gammas2(idx).val, multi_indices);
end
fprintf("Constructing Inverse Map: %d \n",toc)

tic
rebuildSamples = mvnrnd([0,0], eye(2), numSamples);
for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder, idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder, idx);
    else
        multi_indices = genNoCrossMI(polyOrder, idx);
    end
    SAMPLES2(:,idx) = tMAP_vectorized(...
        rebuildSamples(:,1:idx),...
        betas2(idx).val,...
        multi_indices...
        ); 
end

for idx = 1:numDim
    if multi_index_type == "TO"
        multi_indices = genTotalOrderMI(polyOrder, idx);
    elseif multi_index_type == "NM"
        multi_indices = genNoMixedMI(polyOrder, idx);
    else
        multi_indices = genNoCrossMI(polyOrder, idx);
    end
    SAMPLES3(:,idx) = tMAP_vectorized(...
        SAMPLES2(:,1:idx),...
        betas(idx).val,...
        multi_indices...
        ); 
end
fprintf("Computing Inverted Samples: %d \n",toc)
hold on
figure(30)
scatterhist(SAMPLES3(:,1), SAMPLES3(:,2))
title("Low Discrepancy Reconstructed Targets")

%% Plotting More Stuff


figure(4)
[f, xi] = ksdensity(samples);
[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);
hold on
scatter(samples(:,1), samples(:,2), 2.5, 'filled')
hold on
[C,h] = contour(X,Y,Z, 20);
colormap(hot)
h.LineWidth = 2.0;
xlim([1,8])
ylim([1,8])

figure(5)
title("Low-Discrepancy Reconstructed Samples")
[f, xi] = ksdensity(SAMPLES);
[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);
hold on
scatter(SAMPLES(:,1), SAMPLES(:,2), 2.5, 'filled')
hold on
[C,h] = contour(X,Y,Z, 20);
h.LineWidth = 2.0;
xlim([1,8])
ylim([1,8])

figure(6)
uniSamples = ones(5000, 2) + 7*rand(5000, 2);
[f, xi] = ksdensity(uniSamples);
[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);
hold on
scatter(uniSamples(:,1), uniSamples(:,2), 2.5, 'filled')
hold on
[C,h] = contour(X,Y,Z, 20);
colormap(hot)
h.LineWidth = 2.0;
xlim([1,8])
ylim([1,8])

figure(7)
ss = sobolset(2);
uniSamples = ones(5000, 2) + 7*ss(1:5000, :);
[f, xi] = ksdensity(uniSamples);
[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);
hold on
scatter(uniSamples(:,1), uniSamples(:,2), 2.5, 'filled')
hold on
[C,h] = contour(X,Y,Z, 20);
colormap(hot)
h.LineWidth = 2.0;
xlim([1,8])
ylim([1,8])

%% Begin Function Definitions

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

function theta = random_sampler_2(numSamples)
    %theta = zeros([2,numSamples]);
    uniSamples = rand([numSamples,2]);
    refSamples = uniformTOgaussian(uniSamples)';
    r_1 = refSamples(1,:);
    r_2 = refSamples(2,:);
    theta = [1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)] *...
        [r_1; cos(r_1) + 0.5*r_2];
    theta = 5 * ones(size(theta')) + theta';
end