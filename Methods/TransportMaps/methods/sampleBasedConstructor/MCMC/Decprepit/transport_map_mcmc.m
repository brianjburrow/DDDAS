classdef transport_map_mcmc < handle
    %% transport_map_mcmc class
    % Implements Transport Map Accelerated Markov Chain Monte Carlo from
    % Matt Parno's PhD dissertation (Massachusetts Institute of Technology)
    properties
        %% Sample Space Properties
        numDim                                                              % Number of Dimensions
      
        %% MCMC Properties
        adaptGap                                                            % Adapt map every "adaptGap" samples
        samples                                                             % Container for MCMC Samples
        currentIter                                                         % Current chain iteration
        numTotalSamples                                                     % Number of samples desired
        numOriginalSamples                                                  % Number of samples in the initial MCMC chain, used for plotting only
        
        %% TrasportMapProperties
        gammas                                                              % Transport Map Parameters
        betas                                                               % Inverse Transport Map Parameters
        multi_index_type                                                    % Takes values: "TO", "NC", "NM"
        multi_index_func                                                    % Generates Multi-index sets
        multi_indices                                                       % NumDim different multi_index sets
        polynomialOrder                                                     % Max Order of Polynomial Expansion
        sortIndices
        %% Transport Map Objective Func Properties
        F                                                                   % Contains Hermite Polynomial Expansions
        G                                                                   % Contains derivatives of Hermite Polynomial Expansions
        F_mult                                                              % F'*F, stored for efficient updates
        
        %% Transport Map Constraint Func Properties
        dmin                                                                % Small non-zero value.  Hardcoded 10^-8 as in Parno Dissertation
        Kr                                                                  % Small non-zero value.  Hardcoded 10^-4 as in Parno Dissertation
        %% Transport Map Convergence Properties
        
        %% Sampling Distributions
        propSampler                                                         % Generates Proposal Samples
        
        %% Log Probability Density Functions
        log_target_pdf                                                      % log Target of the Markov Chain
        
        %% Optimizers
        optimizer_options                                                   % Matlab optimoptions obj., see template_tm_mcmc.m file
    end
    
    methods
        %% Initialize the transport_map_mcmc class
        function obj = transport_map_mcmc(mcmcProp, tmProp, opt, tarPdf)
           [obj.currentIter, obj.numDim] = size(mcmcProp.samples);         
           obj.adaptGap = mcmcProp.adaptGap;
           obj.numTotalSamples = mcmcProp.numTotalSamples;
           obj.samples = zeros([obj.numTotalSamples, obj.numDim]);
           obj.samples(1:obj.currentIter, :) = mcmcProp.samples;
           obj.numOriginalSamples = obj.currentIter;
           
           obj.multi_index_type = tmProp.multi_index_type;
           obj.polynomialOrder = tmProp.polynomialOrder;
           obj.set_multi_index_func();
           obj.gen_multi_index_sets();
           
           obj.log_target_pdf = tarPdf.log_target_pdf;
           obj.optimizer_options = opt;
           obj.dmin = 10^-8;                                                % Hardcoded value from Matt Parno Dissertation
           obj.Kr = 10^-4;                                                  % Hardcoded value from Matt Parno Dissertation
           fprintf("Initializing Objective Function \n")
           for idx = 1:obj.numDim
               obj.initialize_objective_func_props(idx);
           end
           
           for dmx = 1:obj.numDim      
                value = min(find(obj.multi_indices(dmx).val(:,dmx) > 0));
                MD = length(obj.multi_indices(dmx).val(:,1));
                init = zeros(1, MD);
                init(value) = 1;
                objF = @(xx) obj.objectiveFunc(xx, ...                      % Initialize objective function
                    obj.F_mult(dmx).val, ...
                    obj.G(dmx).val,...
                    init);                                                  % Pass identity map
                const = @(xx) obj.constraintFunc(...                        % Initialize constraint function
                    xx,...
                    obj.G(dmx).val);                                                  
                fprintf("Optimizing Map Component %d \n", dmx)
                obj.gammas(dmx).val = fmincon(...                           % Perform optimization
                                objF, ...
                                init,...                                    % Use previous iteration as initial condition                        
                                [],[],[],[],[],[],...
                                const, ...
                                obj.optimizer_options...
                                );
           end
           fprintf("Computing Reference Samples \n")
           refSamps = zeros([obj.currentIter, obj.numDim]);
            for dmx = 1:obj.numDim
                refSamps(:,dmx) = tMAP_vectorized(...
                    obj.samples(1:obj.currentIter,1:dmx),...
                    obj.gammas(dmx).val,...
                    obj.multi_indices(dmx).val...
                    ); 
            end
            fprintf("Optimizing Inverse Map \n")
            for dmx = 1:obj.numDim
                tarSamples = construct_vandermond(...
                    refSamps,...
                    ones(size(obj.gammas(dmx).val)),...
                    obj.multi_indices(dmx).val);
                obj.betas(dmx).val = tarSamples \ obj.samples(1:obj.currentIter,end);
            end
           fprintf("Beginning Markov Chain Monte Carlo Simulation \n")
           obj.runMCMC();
        end
        %% Functions to Run MCMC
        
        function obj = runMCMC(obj)
            fprintf("Running MCMC \n")
            %% Initialize MCMC run
            numD = obj.numDim;
            numInitSamp = obj.currentIter;
            totalSamp = obj.numTotalSamples;
            samps = obj.samples;
            Ku = obj.adaptGap;
            bet = obj.betas;
            mi = obj.multi_indices;
            propTheta = zeros([1,numD]);
            r_k = zeros([1, obj.numDim]);                                   % Initialize the start of the chain by computing the current reference sample
            for dmx = 1:obj.numDim
                r_k(:,dmx) = tMAP_vectorized(...                            % Map last sample from original MCMC chain to the reference distribution
                    samps(obj.currentIter,1:dmx),...
                    obj.gammas(dmx).val,...
                    obj.multi_indices(dmx).val...
                    ); 
            end
            %% Begin Generating Samples
            for idx = (numInitSamp + 1):totalSamp
                r_p = mvnrnd(zeros([1,numD]), 2*eye(numD));                 % Generate proposal sample from reference distribution
                for dmx = 1:numD
                    propTheta(1,dmx) = tMAP_vectorized(r_p,...              % Map reference proposal to approximate target distribution
                        bet(dmx).val,...
                        mi(dmx).val);
                end
                alpha = obj.logAcceptanceProb(...                           % Compute Acceptance Probability
                    samps(idx-1,:),...
                    propTheta,...
                    r_p,...
                    r_k);
                if rand <= exp(alpha)                                       % Take exp b/c we have log(alpha);
                    samps(idx,:) = propTheta;                               % Accepted Proposal, store current iteration
                    r_k = r_p;                                              % Update the reference proposal;
                else
                    samps(idx,:) = samps(idx-1, :);                         % Rejected Proposal, store previous iteration
                end
                
                if mod(idx, Ku) == 0                                        % Update the transport map
                    fprintf(...
                        "Updating Transport Map at Iteration %d \n",...
                        idx);
                    obj.currentIter = idx;
                    for dmx = 1:numD
                        value = min(find(obj.multi_indices(dmx).val(:,dmx) > 0));
                        MD = length(obj.multi_indices(dmx).val(:,1));
                        init = zeros(1, MD);
                        init(value) = 1;
                        obj.update_objective_func_props(dmx);
                        objF = @(xx) obj.objectiveFunc(xx, ...              % Initialize objective function w/updated F_mult, G
                            obj.F_mult(dmx).val, ...
                            obj.G(dmx).val,...
                            init);
                        const = @(xx) obj.constraintFunc(...                % Initialize constraint function
                            xx,...
                            obj.G(dmx).val...
                            );
                        obj.gammas(dmx).val = fmincon(...                   % Perform optimization
                                        objF, ...
                                        obj.gammas(dmx).val,...             % Can substitute obj.gammas(dmx).val to use previous iteration as initial condition                        
                                        [],[],[],[],[],[],...
                                        const, ...
                                        obj.optimizer_options...
                                        );
                    end
                    fprintf("Recomputing Reference Samples \n")
                    refSamps = zeros([idx, obj.numDim]);                    % Recompute all reference samples
                    for dmx = 1:obj.numDim
                        refSamps(:,dmx) = tMAP_vectorized(...
                            samps(1:idx,1:dmx),...
                            obj.gammas(dmx).val,...
                            obj.multi_indices(dmx).val...
                            ); 
                    end
                    fprintf("Inverting Map \n")
                    for dmx = 1:numD                                        % Update Inverse Map with new reference samples
                        tarSamples = construct_vandermond(...
                            refSamps,...
                            ones(size(obj.gammas(dmx).val)),...
                            obj.multi_indices(dmx).val);
                        obj.betas(dmx).val = tarSamples \ samps(1:idx,1:dmx);
                    end
                end
            end
            obj.currentIter = idx;
            obj.samples = samps;
        end
        
        %% Function Handle Targets
        % Base functions that we alter to create function handles
        % with certain inputs constant
        function [value,value_der] = objectiveFunc(obj, params, f_mult, g,...
                identityMap)
            params = params';                                               % Params is either gammas or betas.  Depends on what is passed.

            %% Compute Objective Value
            preF = params'*f_mult;
            G_mult = g*params;
            regDiff = params - identityMap;
            regularizer = norm(regDiff)^2;
            value = 0.5*preF*params - ...
                ones([1,obj.currentIter]) * log(G_mult) +...
                obj.Kr*regularizer;

            %% Compute Derivative
            if nargout > 1
                value_der = zeros(size(params));
                postF = f_mult*params;
                for idx = 1:length(params')                                 % Vectorize this if possible
                    value_der(idx) = 0.5*(postF(idx) +...
                        preF(idx)) -...
                        sum(g(:,idx)./G_mult) + ...
                        2*regDiff(idx);
                end
            end
            
        end
        
        function [value, eqCon, DC, DCeq] = constraintFunc(obj, gamma, g)
            %% Variables
            % Samples:       from target distribution
            % Gamma:         parameters of the transport map
            % Index:         Take derivatives w.r.t
            % Multi_indices: Determines polynomial expansion
            eqCon = [];
            %% Construct F,G Matrices
            value = obj.dmin*ones([obj.currentIter,1]) - g*gamma';
            if nargout > 2
                DCeq = [];
                DC = -g';
            end
        end
        
        %% Utility Functions
        function log_p = log_prop_pdf(~, input, condition)
            log_p = logmvnpdf(input, condition, eye(length(input)));
        end
        
        function obj = set_multi_index_func(obj)
            type = obj.multi_index_type;
            if type == "TO"
                obj.multi_index_func = @(po, mi) genTotalOrderMI(po, mi);
            elseif type == "NM"
                obj.multi_index_func = @(po, mi) genNoMixedMI(po, mi);
            elseif type == "NC"
                obj.multi_index_func = @(po, mi) genNoCrossMI(po, mi);
            else
                fprintf("No multi-index function for type %s", type);
            end
        end
        
        function obj = gen_multi_index_sets(obj)
            for idx = 1:obj.numDim
                obj.multi_indices(idx).val = obj.multi_index_func(...
                    obj.polynomialOrder, ...
                    idx);
            end
        end
        
        function obj = initialize_objective_func_props(obj, dimension)
            %% Construct F,G Matrices
            samps = obj.samples(1:obj.currentIter, 1:dimension);
            MI = obj.multi_indices(dimension).val;
            K = length(samps(1:obj.currentIter,1));
            M = length(MI(:,1));
            f    = zeros([K, M]);
            g    = zeros([K, M]);

            mvp  = @(xx, yy) multivariatePolynomial_vectorized(xx, yy);
            pmvp = @(xx, yy) partialMVP_vectorized(xx, yy, dimension);

            for dmx = 1:M
                f(:, dmx) = mvp(...
                                    samps,...
                                    MI(dmx, 1:dimension)...
                                    );

                g(:, dmx) = pmvp(...
                                    samps,...
                                    MI(dmx,1:dimension)...
                                    );
            end
            obj.F(dimension).val = f;
            obj.G(dimension).val = g;
            
            obj.F_mult(dimension).val = obj.F(dimension).val'*obj.F(dimension).val;
        end
        
        function obj = update_objective_func_props(obj, dimension)
            %% Update F, G, F_mult
            prevIter = length(obj.F(dimension).val(:,1)) + 1;
            
            samps = obj.samples(prevIter:obj.currentIter, 1:dimension);
            MI = obj.multi_indices(dimension).val;
            K = length(samps(:,1));
            M = length(MI(:,1));
            
            f    = zeros([K, M]);
            g    = zeros([K, M]);
            
            mvp  = @(xx, yy) multivariatePolynomial_vectorized(xx, yy);
            pmvp = @(xx, yy) partialMVP_vectorized(xx, yy, dimension);

            for dmx = 1:M
                f(:, dmx) = mvp(...
                                    samps,...
                                    MI(dmx, 1:dimension)...
                                    );

                g(:, dmx) = pmvp(...
                                    samps,...
                                    MI(dmx,1:dimension)...
                                    );
            end
            obj.F(dimension).val = [obj.F(dimension).val; f];
            obj.G(dimension).val = [obj.G(dimension).val; g];
            obj.F_mult(dimension).val = obj.F_mult(dimension).val + f'*f;  
        end
        
        function alpha = logAcceptanceProb(obj, currentTheta, propTheta, r_p, r_k)
            logP = obj.log_target_pdf(propTheta) +...
                obj.log_prop_pdf(r_k, r_p) - ...
                obj.log_target_pdf(currentTheta) - ...
                obj.log_prop_pdf(r_p, r_k); 
            alpha = min(log(1), logP);
        end     
        
        function plotMarginalHists(obj)
            %% Generate Samples from Transport Map to show Proposal Samples
            randSamples = normrnd(0, 1, [10000, 1]);
            targetSamples = zeros(size(randSamples));
            for dmx = 1
                targetSamples(:,dmx) = tMAP_vectorized(...
                                            randSamples(:,1:dmx),...
                                            obj.betas(dmx).val,...
                                            obj.multi_indices(dmx).val...
                                            );                      
            end
                
            n = get(gcf,'Number');
            for idx = 1:obj.numDim
                figure(n + idx)
                subplot(2,2,1)
                histogram(obj.samples(1:obj.numOriginalSamples, idx))
                title("Original MCMC Samples")

                subplot(2,2,2)
                histogram(obj.samples(:, idx))
                title("Total MCMC Samples")

                subplot(2,2,[3,4])
                histogram(targetSamples(:,idx))
                title("Transport Map Generated Samples")
            end
        end
        function plotTrace(obj)
            n = get(gcf,'Number');
            figure(n + 1)
            for idx = 1:obj.numDim
                subplot(obj.numDim, 1, idx)
                scatter(1:obj.currentIter, obj.samples(:, idx))
                tit = sprintf("Trace of MCMC Samples: Dimension %d", idx);
                title(tit)
            end
        end
        
        function plotReferenceDistributions(obj)
            n = get(gcf, 'Number');
            fprintf("Computing Reference Samples \n")
            refSamps = zeros([obj.currentIter, obj.numDim]);                % Recompute all reference samples
            for dmx = 1:obj.numDim
                refSamps(:,dmx) = tMAP_vectorized(...
                    obj.samples(1:obj.currentIter,1:dmx),...
                    obj.gammas(dmx).val,...
                    obj.multi_indices(dmx).val...
                    ); 
            end
            randNormal = normrnd(0, 1, [obj.currentIter, obj.numDim]);
            reconstructedSamples = zeros([obj.currentIter, obj.numDim]);
            for dmx = 1:obj.numDim
                reconstructedSamples(:,dmx) = tMAP_vectorized(...
                    randNormal(:,1:dmx),...
                    obj.betas(dmx).val,...
                    obj.multi_indices(dmx).val...
                    ); 
            end
            for dmx = 1:obj.numDim
                figure(n + 1)
                subplot(3,1,1)
                histogram(obj.samples(:,dmx))
                title("MCMC Samples")
                subplot(3,1,2)
                histogram(refSamps(:,dmx))
                title("MCMC Samples mapped to Reference")
                subplot(3,1,3)
                histogram(reconstructedSamples(:,dmx))
                title("Transport Map Approximation")
            end
        end
    end
end