function input = L2oFilter_up(input, truncate, numStartSamples, itr, knnTrue)
    % Current best: on Nov 1, 2018
    %% Input
    % L2O.samples = propSamples';                                           % Offline library of state vectors (nDimension X nSamples)
    %L2O.data = cappe_benchmark.GetDataC(0);                                % Observation from the system at the current time step
    %L2O.stateEvolution = @(xx, yy)                                         % State evolution function, first input is the state, second input are operating conditions
    %L2O.oper = 0;                                                          % Operating conditions.  
    %L2O.output = output';                                                  % forward model evaluations of the Offline library of state vectors
    %L2O.R = measurementNoise;                                              % Measurement Noise covariance matrix
    %L2O.estimate = 0;                                                      % Previous estimate of the state (or initial estimate of the state)
    %L2O.err = processNoise;                                                % Previous estimate of the error covariance matrix
    %% Input sizing
    %
    %% Unpackage Input object
    samples = input.samples;
    data = input.data;
    stateEvolution = input.stateEvolution;
    oper = input.oper;
    output = input.output;
    R = input.R;
    
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
    %% Calc array sizes
    numSamples = length(samples(1,:));
 
    %% Determine Weights
    mins = min([samples, newSamples], [], 2);                               % Find lower coordinates of a bounding box
    maxs = max([samples, newSamples], [], 2);                               % Find upper coordinates of a bounding box
    normedProps = (samples - mins)./(maxs - mins);                          % Use bounding box to transform problem to [0,1]^d hypercube
    normedTargets = (newSamples - mins)./(maxs - mins);                     % Use bounding box to transform problem to [0,1]^d hypercube
    [normedProps, propIdx] = sortrows(normedProps');                        % Sort the proposal samples
    [normedTargets, ~] = sortrows(normedTargets');                          % Sort the target samples
    output = output(:,propIdx);                                             % Sort the offline library forward model evaluations to match the new proposal ordering
    samples = samples(:,propIdx);                                           % Sort the unnormalized proposal samples
    
    w = ones(numSamples,1)./numSamples;                                     % Initialize a feasible set of weights
    H = matvecProd(normedProps, w);                                         % Compute the hadamard product matrix
    B = buildXYb(normedProps, normedTargets);                               % Compute the B vector
    feval = zeros(itr+1, 1);                                                % Initialize the objective function eval vector
    [weights, ~] = frankwolfe(normedProps, B, w, H, feval, itr);            % Compute the set of L2-optimal importance weights
     %EFF = sum(weights)^2 / (sum(weights.^2));                             % Compute the effective sample size if desired
     %fprintf('Effective Sample Size %0.5f \n', EFF);
    
    %% Predict the mean and variance of the state prior
    mu = mean(newSamples,2);                                                % Compute the mean of the target samples
    covariance = cov(newSamples);                                           % Compute the covariance of the target samples
    
    %% Predict the observation
    obs = zeros(length(data(:,1)), 1);
    for idx = 1:numSamples
        obs = obs + weights(idx) * output(:,idx);                           % Predict the observation using the weighted forward model evaluations
    end
    
    %% Compute the innovation covariance
    P_yy = R;
    for idx = 1:numSamples
        P_yy = P_yy + ...
            weights(idx) * (output(:,idx) - obs)*(output(:,idx) - obs)';
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
    input.estimate = x;                                                     % New estimate of the state
    input.err = diag(P);                                                    % New estimate of the uncertainty
    
end