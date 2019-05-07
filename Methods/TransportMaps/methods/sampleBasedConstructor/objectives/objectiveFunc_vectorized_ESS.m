function [value,value_dir] = objectiveFunc_vectorized_ESS(gamma, samples, multi_indices, d,...
    ESS_func, pointSet, polyOrder)
    %% Variables
    % Samples:       from target distribution
    % Gamma:         parameters of the transport map
    % Multi_indices: Determines polynomial expansion
    K = length(samples(:,1));
    M = length(multi_indices(:,1));
    gamma = gamma';
    %% Construct F,G Matrices
    F    = zeros([K, M]);
    G    = zeros([K, M]);

    mvp  = @(xx, yy) multivariatePolynomial_vectorized(xx, yy);
    pmvp = @(xx, yy) partialMVP_vectorized(xx, yy, d);

    for dmx = 1:M
        F(:, dmx) = mvp(...
                            samples(:,1:d),...
                            multi_indices(dmx, :)...
                            );

        G(:, dmx) = pmvp(...
                            samples(:,1:d),...
                            multi_indices(dmx,:)...
                            );
    end
    
    %% Compute Inverse Map
    refSamples = zeros(size(samples));

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
    
    refSamples = pointSet;
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
    
    worstCaseESS = ESS_func(SAMPLES);
    %% Compute Objective
    F_mult = F'*F;
    preF = gamma'*F_mult;
    G_mult = G*gamma;
    value = 0.5*preF*gamma - ones([1,K]) * log(G_mult) - worstCaseESS;

    %% Compute Derivative
    if nargout > 1000
        value_dir = zeros(size(gamma));
        postF = F_mult*gamma;
        for idx = 1:length(gamma')
            value_dir(idx) = 0.5*(postF(idx) + preF(idx)) - sum(G(:,idx)./G_mult);
%             for dmx = 1:K
%                 value_dir(idx) = value_dir(idx) - G(dmx,idx)/(G_mult(dmx,:));
%             end
        end
    end
end