function weights = plot_l2o_update(proposal_samples, target_samples)
    %% Proposal samples
    samples = proposal_samples;
    %% Target Samples
    newSamples = target_samples;
    %% Termination Criterion
    itr = 10000;

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
    
     subplot(3,3,9)
     stairs(samples', swCDF, 'black')

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
end