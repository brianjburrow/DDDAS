function [value,value_dir] = objectiveFunc_vectorized(gamma, samples, multi_indices, d)
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
    %% Compute Objective
    F_mult = F'*F;
    preF = gamma'*F_mult;
    G_mult = G*gamma;
    value = 0.5*preF*gamma - ones([1,K]) * log(G_mult);

    %% Compute Derivative
    if nargout > 1
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