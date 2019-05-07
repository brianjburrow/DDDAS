classdef ActiveSubspace < handle
    %% Implements Active Subspaces
    % Author: Brian Burrows - Brianjburrow@gmail.com
    % Projects an input space onto a lower dimension.
    % That lower dimension is found by looking at the average squared 
    % gradient, where the average is taken with respect to a particular
    % probability distribution.
    
    % INPUT:                                                                type:             Description
    % gradFunc:                                                             function handle:  takes in a single input of nDimensions
    % nSamples:                                                             Scalar:           the number of Monte Carlo samples used to compute the active subspace (the higher the more accurate)
    % randomNumGenerator:                                                   function handle:  takes in a single input, the number of samples desired.  Returns samples in the form nSamples x nDimensions 
    % eigCutoff:                                                            Scalar:           Used to choose the active subspace dimensionality.  Range [0,1].
    % nDimensions:                                                          scalar:           Dimensionality of input space
    
    % Parameters (proprties):                                               type:             Description
    properties
        % Input properties
        gradFunc                                                            % See above
        nSamples
        randomNumGenerator
        eigCutoff
        nDimensions
        % Computed Properties
        activeDim                                                           % Scalar:         Number of dimensions in the active subspace
        C                                                                   % Matrix:         nDimensions x nDimensions describing the average derivative functional
        gradients                                                           % Matrix:         nDimensions x nSamples: contains nSample evaluations of gradient function
        Ds                                                                  % Matrix:         nDimensions x nDimensions: contains diagonal matrix of eigenvalues from eigendecomposition
        Vs                                                                  % Matrix:         nDimensions x nDimensions: contains right eigenvectors from eigendecomosition
        W1                                                                  % Matrix:         nDimensions x activeDim:   contains active subspace basis vectors
        W2                                                                  % Matrix:         nDimensions x nDimensions - activeDim: contains inactive subspace basis vectors
        DsNormed                                                            % Matrix:         size(Ds) : contains normalized eigenvalues for computing active subspace dimension
    end
    
    methods
        function obj = ActiveSubspace(gradFunc, nSamples,...
                nDimensions, randomNumGenerator,...
                eigCutoff)
            %% Initialize the ActiveSubspace object
            obj.gradFunc = gradFunc;
            obj.nSamples = nSamples;
            obj.randomNumGenerator = randomNumGenerator;
            obj.nDimensions = nDimensions;
            obj.eigCutoff = eigCutoff;
            obj.ComputeGradients();
            obj.ComputeC();
            obj.ComputeActiveSubspace();
        end
        function obj = ComputeGradients(obj)
            %% Precompute gradient matrix used for Monte carlo approx of obj.C
            obj.gradients = zeros([obj.nDimensions, obj.nSamples]);
            samples = obj.randomNumGenerator(obj.nSamples);
            obj.gradients = obj.gradFunc(samples);
        end
        function obj = ComputeC(obj)
            %% Compute average gradient functional using Monte Carlo Integration
            obj.C = zeros(obj.nDimensions);
            for iSample = 1:obj.nSamples
                obj.C = obj.C + obj.gradients(:, iSample)*...
                    obj.gradients(:, iSample)';
            end
            obj.C = obj.C / obj.nSamples;
        end
        function obj = ComputeActiveSubspace(obj)
            %% Find the active subspace dimension and related basis vectors
            [V, D, ~] = eig(obj.C);
            disp('V')
            disp(V)
            [~, ind] = sort(diag(D), 'descend');
            obj.Ds = D(ind, ind);
            obj.Vs = V(:, ind);
            obj.DsNormed = obj.Ds ./ sum(diag(obj.Ds));
            cumDs = cumsum(diag(obj.DsNormed));
            obj.activeDim = find(cumDs > obj.eigCutoff, 1, 'first');
            obj.W1 = obj.Vs(:,1:obj.activeDim);
            obj.W2 = obj.Vs(:, obj.activeDim + 1:end);
        end
        
        function y = ProjectSamples(obj, X)
            %% Project samples down into lower dimensional space
            % This function is intended to be used after computing the 
            % active subspace.  
            y = (obj.W1' * X')';
        end
        
        function x = invertProjection(obj, Y)
            x = ((obj.W1')^(1) * Y')';
        end
    end
    
end