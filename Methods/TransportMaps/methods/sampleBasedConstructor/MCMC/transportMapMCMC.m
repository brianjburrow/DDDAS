classdef transportMapMCMC < handle
    properties
        samples
        polyOrder
        gammas
        adaptSpacing
        numTotalSamples
        currentIter
        log_likelihood
        log_prior
        numInitSamples
        numDim
        betas
        log_map_derivative
        log_target_density
        chainInit
        propVar
        multi_index_type
    end
    
    methods
        function obj = transportMapMCMC(log_likelihood, log_prior,...
                initSample, polyOrder, adaptSpacing, numTotalSamples...
                )
            disp('in')
            obj.numInitSamples = length(initSample(:,1));
            obj.samples = zeros([numTotalSamples, length(initSample(1,:))]);
            obj.log_likelihood = log_likelihood;
            obj.log_prior = log_prior;
            obj.samples(1:obj.numInitSamples,:) = initSample;
            obj.polyOrder = polyOrder;
            obj.adaptSpacing = adaptSpacing;
            obj.numTotalSamples = numTotalSamples;
            obj.numDim = length(initSample(1,:));
            
            obj.log_map_derivative = @(xx, yy, zz, tmz) log(tMAP_partial_vectorized(...
                xx, yy, zz, tmz...
                ));
            
            obj.log_target_density = @(xx) log_likelihood(xx) + log_prior(xx); % xx is the variable of interest, data should be hard coded into likelihood
            obj.multi_index_type = "TO";
            obj.initializeMap();
            obj.runMCMC();
        end
        
        function obj = initializeMap(obj)
            samps = obj.samples(1:obj.numInitSamples,:);  
            polyOrd = obj.polyOrder;
            mi_type = obj.multi_index_type;
            for dmx = 1:obj.numDim
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  dmx);
                elseif mi_type == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  dmx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  dmx);
                end
                value = min(find(multi_indices(:,dmx) > 0));
                MD = length(multi_indices(:,1));
                init = zeros(1, MD);
                init(value) = 1;
                if dmx > 1
                    numGamPrev = length(obj.gammas(dmx-1).val);
                    init(1:numGamPrev) = obj.gammas(dmx-1).val;
                end
                obj.gammas(dmx).val = optimizeMap(samps(:,1:dmx),...
                                                 [],...
                                                 multi_indices, init);
            end
            obj.initializeInverse(samps);
        end
        
        function obj = updateMap(obj, currentNumSamples)
            samps = obj.samples(1:currentNumSamples,:);  
            polyOrd = obj.polyOrder;
            mi_type = obj.multi_index_type;
            for dmx = 1:obj.numDim
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  dmx);
                elseif mi_type == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  dmx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  dmx);
                end
                value = min(find(multi_indices(:,dmx) > 0));
                MD = length(multi_indices(:,1));
                init = zeros(1, MD);
                init(value) = 1;
                if dmx > 1
                    numGamPrev = length(obj.gammas(dmx-1).val);
                    init(1:numGamPrev) = obj.gammas(dmx-1).val;
                end
                obj.gammas(dmx).val = optimizeMap(...
                                                 samps(:,1:dmx),...
                                                 [],...
                                                 multi_indices, ...
                                                 obj.gammas(dmx).val...
                                                 );
            end
            obj.initializeInverse(samps);
        end
        
        function obj = initializeInverse(obj, samps)
            %% Determine Reference Samples
            polyOrd = obj.polyOrder;
            gam = obj.gammas;
            refSamples = zeros(size(samps));
            mi_type = obj.multi_index_type;
            for idx = 1:obj.numDim
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  idx);
                elseif mi_type == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  idx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  idx);
                end
                refSamples(:,idx) = tMAP_vectorized(samps(:,1:idx), gam(idx).val, multi_indices); 
            end
            %% Determine Inverse Map Parameters
            for idx = 1:obj.numDim
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  idx);
                elseif mi == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  idx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  idx);
                end
                
                obj.betas(idx).val = invertMap(samps(:,1:idx), refSamples(:,1:idx), gam(idx).val, multi_indices);
            end
        end
        
        function obj = runMCMC(obj)
            numD = obj.numDim;
            numInitSamp = obj.numInitSamples;
            totalSamp = obj.numTotalSamples;
            samps = obj.samples;
            Ku = obj.adaptSpacing;
            logAcceptProb = @(xx, yy) obj.logAcceptanceProb(xx, yy);       % put in local scope for speed
            for idx = (numInitSamp + 1):totalSamp
                prop = mvnrnd(zeros([1,numD]), 2*eye(numD));
                alpha = logAcceptProb(samps(idx-1,:), prop);
                if rand <= alpha
                    samps(idx,:) = prop;
                else
                    samps(idx,:) = samps(idx-1, :);
                end
                
                if mod(idx, Ku) == 0
                    fprintf("Updating Transport Map at Iteration %d", idx)
                    obj.updateMap(idx);
                end
            end
            obj.samples = samps;
        end
        
        function log_p = log_target(obj, input, refSample)
            log_lik = obj.log_target_density;
            D = obj.numDim;
            log_p = log_lik(input);
            bet = obj.betas;
            polyOrd = obj.polyOrder;
            for idx = 1:D
                mi_type = obj.multi_index_type;
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  idx);
                elseif mi_type == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  idx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  idx);
                end
                log_p = log_p + tMAP_partial_vectorized(refSample, bet(idx).val, multi_indices, idx);
            end
        end
        
        function log_p = log_propPdf(~, input, condition)
            log_p = logmvnpdf(input, condition, eye(length(input)));
        end
        
        function alpha = logAcceptanceProb(obj, currentTheta, r_p)
            
            %% Determine Reference Samples
            polyOrd = obj.polyOrder;
            gam = obj.gammas;
            bet = obj.betas;
            mi_type = obj.multi_index_type;
            r_k = zeros(size(currentTheta));
            propTheta = zeros(size(currentTheta));
            for idx = 1:obj.numDim
                if mi_type == "TO"
                    multi_indices = genTotalOrderMI(polyOrd,  idx);
                elseif mi_type == "NM"
                    multi_indices = genNoMixedMI(polyOrd,  idx);
                else
                    multi_indices = genNoCrossMI(polyOrd,  idx);
                end
                r_k(:,idx) = tMAP_vectorized(currentTheta(:,1:idx), gam(idx).val, multi_indices); 
                propTheta(:,idx) = tMAP_vectorized(r_p(:,1:idx), bet(idx).val, multi_indices);
            end
            
            logP = obj.log_target(propTheta, r_p) +...
                obj.log_propPdf(r_k, r_p) - ...
                obj.log_target(currentTheta, r_k) - ...
                obj.log_propPdf(r_p, r_k);
            
            alpha = min(log(1), logP);
        end       
    end
end