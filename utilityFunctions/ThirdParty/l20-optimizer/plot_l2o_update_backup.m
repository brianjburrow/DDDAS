clear all
close all

%% Proposal samples
samples = normrnd(0, 4, [1,25]);
%% Target Samples
newSamples = normrnd(1, 0.5, [1, 20]);
%% Termination Criterion
itr = 200000;

%% Determine Weights
numSamples = length(samples(1,:));
concatSamps = [samples, newSamples];
mins = min([samples, newSamples], [], 2);
maxs = max([samples, newSamples], [], 2);
normedProps = (samples - mins)./(maxs - mins);
normedTargets = (newSamples - mins)./(maxs - mins);
[normedProps, propIdx] = sortrows(normedProps');
[normedTargets, ~] = sortrows(normedTargets');
samples = samples(:,propIdx);

w = ones(numSamples,1)./numSamples;
H = matvecProd(normedProps, w);
B = buildXYb(normedProps, normedTargets);
% Construct B vector
feval = zeros(itr+1, 1);
[weights, feval] = frankwolfe(normedProps, B, w, H, feval, itr);
 EFF = sum(weights)^2 / (sum(weights.^2));
 
 
 swCDF(1) = weights(1);
 for idx = 2:length(weights)
     swCDF(idx) = swCDF(idx - 1) + weights(idx);
 end
 
 figure(1)
 xlabel('Sample Space')
 hold on
 ylabel('P(x < X)')
 ecdf(samples);
 ecdf(newSamples);
 stairs(samples', swCDF, 'black')
 ylim([0,1])
 legend("Proposal", "Target", "Weights")
 
 function prob = evalCDF(value, cdfArr, cdfArr2)
    for tmz = 1:length(value)
        iiddxx = max(cdfArr(cdfArr2' <= value(tmz)));
        if isempty(iiddxx)
            prob(tmz) =  0;
        else
            prob(tmz) =  iiddxx;
        end
    end
end