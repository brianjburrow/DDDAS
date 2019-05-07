function input = BootstrapFilter(input)
    %% Particle Filter(input)
    % Input
    % input.x = state vector (column vector N x 1)
    
    %% Unpack Input Object
    particles = input.particles;                                            % Set of states, nSamples X nDimensions
    weights = input.weights;                                                % Corresponding set of weights for each particle
    resampling = input.resampling;                                          % Binary (0 or 1) that decides whether or not to resample for a particular iteration
    likelihood = input.likelihood;                                          % likelihood (probability density function)
    proposalSampler = input.proposalSampler;                                % a random number generator for the proposal distribution
    oper = input.oper;                                                      % operating parameters (if needed)
    data = input.data;                                                      % measurement from the system at the current time step

    
    N = length(particles(1,:));                                             % Compute the number of particles we have
    %% Resampling
    if resampling        
        disp("RESAMPLING")
        tempPart = zeros(size(particles));
        count = 1;
        r = mnrnd(N, weights);                                              % Multinomial resampling strategy (Not recommended in general)
                                                                            % May be useful to replace this with some other strategy
        for idx = 1:N
            for dmx = 1:r(idx)
                tempPart(:, count) = particles(:,idx);
                count = count + 1;
            end
        end
        weights = (1/N) * ones(size(weights));
        particles = tempPart;
    end
    
    %% Propagate Particles
    rng = @(xx) proposalSampler.random(xx, data, oper);                     % Set the random number generator with a fixed operating condition
    particlesNew = rng(particles);
    likHood = @(xx) likelihood(data, xx);                                   % fix the data in the likelihood function
    likHoodEval = likHood(particlesNew);
    weights = weights.*likHoodEval;
    weights = weights./sum(weights);
    input.ESS = sum(weights)^2 / (sum(weights.^2));                         % Compute the ESS to decide whether to resample
    input.weights = weights;                                                % Update the weights
    input.particles = particlesNew;                                         % Update the particles
    temp = weights * particlesNew';
    input.estimate = temp;                                                  % Update the estimated state
    temp = (particlesNew - temp).^2;
    variance = weights*temp';

    input.err = variance;                                                   % Update the estimated state error
end
    
    
    