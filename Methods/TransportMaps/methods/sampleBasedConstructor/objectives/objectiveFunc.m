function value = objectiveFunc(gamma, samples, multi_indices, d)
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

mvp  = @(xx, yy) multivariatePolynomial(xx, yy);
pmvp = @(xx, yy) partialMVP(xx, yy, d);

for idx = 1:K
    for dmx = 1:M
        F(idx, dmx) = mvp(...
                            samples(idx,1:d),...
                            multi_indices(dmx, 1:d)...
                            );
        
        G(idx, dmx) = pmvp(...
                            samples(idx,1:d),...
                            multi_indices(dmx,1:d)...
                            );
    end
end

%% Compute Objective
value = 0;
value = value + ...
    0.5*gamma'*(F'*F)*gamma - ones([1,K]) * log(G*gamma);