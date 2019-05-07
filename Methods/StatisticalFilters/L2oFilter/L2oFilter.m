function input = L2oFilter(input, truncate)
    % PREFERRED FILE IS L2oFilter_up, which is an updated version of this file with better commenting and cleaner code
    % This is used in a few demo files though, so I'm leaving it rather than recoding everything
    % Current best: on Nov 1, 2018
    %% Input
    % L2O.samples = propSamples';                                           % Offline library of state vectors (nDimension X nSamples)
    %L2O.data = cappe_benchmark.GetDataC(0);                                % 
    %L2O.stateEvolution = @(xx, yy)
    %State 
    %L2O.oper = 0;
    %L2O.output = output';
    %L2O.dataFunc = @(xx, yy) cappe_benchmark.GetDataC(xx);
    %L2O.R = measurementNoise;
    %L2O.estimate = 0;
    %L2O.err = processNoise;                                             
    %% Unpackage Input object
    samples = input.samples;
    data = input.data;
    stateEvolution = input.stateEvolution;
    oper = input.oper;
    output = input.output;
    R = input.R;
    
    itr = 5000;                                                             % Number of Optimization Iterations
    numStartSamples = 2000;
    knnTrue = 1;

    numDimensions = length(samples(:,1));
    num = truncate;                                                         % Decide whether to truncate the distribution
    if num == 0
        startSamples = mvnrnd(input.estimate', diag(input.err), numStartSamples)';
    else
        startSamples = zeros([numDimensions, numStartSamples]);
        for idx = 1:numDimensions
            pd = makedist('Normal');
            pd.mu = input.estimate(idx,1);
            pd.sigma = input.err(idx)^0.5;
            if idx == 1
                t = truncate(pd, 0, 1.0);
            elseif idx == 2
                t = truncate(pd, 0, 1.0);
            else
                t = truncate(pd, 0, 1.0);
            end
            startSamples(idx,:) = random(t,[1,numStartSamples]);
        end
    end
    
    %% Propagate Samples
    newSamples = zeros([numDimensions, numStartSamples]);
    for idx = 1:numStartSamples
        newSamples(:,idx) = stateEvolution(startSamples(:,idx), oper);
    end

     if knnTrue
         [indx, ~] = knnsearch(samples', newSamples','K', 1);
         indx = sort(indx);
         samples = samples(:,indx); 
         output = output(:,indx);
     end
     disp('sizeSamples')
     disp(size(samples));
    %% Calc array sizes
    numSamples = length(samples(1,:));
 
    %% Determine Weights
    concatSamps = [samples, newSamples];
    mins = min([samples, newSamples], [], 2);
    maxs = max([samples, newSamples], [], 2);
    normedProps = (samples - mins)./(maxs - mins);
    normedTargets = (newSamples - mins)./(maxs - mins);
    [normedProps, propIdx] = sortrows(normedProps');
    [normedTargets, ~] = sortrows(normedTargets');
    output = output(:,propIdx);
    samples = samples(:,propIdx);
    
    w = ones(numSamples,1)./numSamples;
    H = matvecProd(normedProps, w);
    B = buildXYb(normedProps, normedTargets);
    % Construct B vector
    feval = zeros(itr+1, 1);
    [weights, feval] = frankwolfe(normedProps, B, w, H, feval, itr);
     EFF = sum(weights)^2 / (sum(weights.^2));
     %fprintf('Effective Sample Size %0.5f \n', EFF);
    
    %% Predict The Mean  Value
% 
%     mu = zeros(numDimensions, 1);
%     for idx = 1:numSamples
%         mu = mu + weights(idx)*samples(:,idx);
%     end
    mu = mean(newSamples,2);
    
    covariance = zeros(length(samples(:,1)));
    for idx = 1:numSamples
        covariance = covariance + weights(idx)*(samples(:,idx) - mu)*(samples(:,idx) - mu)';
    end
    
    %% Predict the observation
    obs = zeros(length(data(:,1)), 1);
    for idx = 1:numSamples
        obs = obs + weights(idx) * output(:,idx);
    end
    
    %% Predict the innovation covariance
    P_yy = R;
    for idx = 1:numSamples
        P_yy = P_yy + weights(idx) * (output(:,idx) - obs)*(output(:,idx) - obs)';
    end
    %% Predict the cross correlation matrix
    for idx = 1:numSamples
        if idx == 1
            P_xy = weights(idx)*(samples(:,idx) - mu)*(output(:,idx) - obs)';
        else
            P_xy = P_xy + weights(idx)*(samples(:,idx) - mu)*(output(:,idx) - obs)';
        end
    end
    
    %% Calculate Kalman Gain Matrix
    K = P_xy * P_yy^-1;
    x = mu + K*(data - obs);
    P = covariance - K*P_yy*K';
    
    %% Repackage Object
    input.estimate = x;
    input.err = diag(P);
    
    
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