function obj_dir = objective_derivative(gamma, samples, multi_indices, d)

%% Variables
% Samples:       from target distribution
% Gamma:         parameters of the transport map
% Multi_indices: Determines polynomial expansion
[K, dim] = size(samples);
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
                        multi_indices(dmx, 1:d)...
                        );

    G(:, dmx) = pmvp(...
                        samples(:,1:d),...
                        multi_indices(dmx,1:d)...
                        );
end

%% Compute Objective
value = zeros([dim, 1]);
F_mult = F'*F;
preF = gamma'*F_mult;
postF = F_mult*gamma;
denom = G*gamma;
left = ones([1,K]) * (1/(denom));
for idx = 1:dim
    vec = zeros(size(gamma));
    vec(idx) = 1;
    value(idx, 1) = 0.5*(preF(idx) + postF(idx)) - left*log(G*vec);
end