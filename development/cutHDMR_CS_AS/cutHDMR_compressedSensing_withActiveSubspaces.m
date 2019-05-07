classdef cutHDMR_compressedSensing_withActiveSubspaces < handle
    %% Constructs a cutHDMR_compressedSensing 
    % Implements compressed sensing based cutHDMR developed by Kaiyu Li and
    % Douglas Allaire in "A Compressed Sensing Approach to Uncertainty Propagation
    % for approximately additive functions"
    
    % Author: Brian Burrows - Brianjburrow@gmail.com
    % Input:                                                                Type :           Description
    % fullModel:                                                            function handle: the model you want to approximate, takes single matrix valued input;
    % order:                                                                scalar:          0, 1, or 2 - controls HDMR expansion order
    % cutLine                                                               vector:          nDimensions x 1
    % nSamples                                                              Scalar: 
    properties
        cutLine                                                             
        cutOutput
        inputSet                                                            
        nSamples
        nDimensions
        targetFunction
        sampleSet
        pointGrid
        polynomialBasisType                                                 %String : "legendre", "hermite", "bessel"
        polyFunc
        inputSampler
        nOrder
        L_struct     
        Psi_struct
        mainEffectSensivitivies
    end
    
    methods
        function obj = cutHDMR_compressedSensing_withActiveSubspaces(fullModel, cutLine,...
                nSamples, nOrder, polyType, inputSampler)
            [~, obj.nDimensions] = size(cutLine);
            obj.nSamples         = nSamples;
            obj.cutLine          = cutLine;
            obj.cutOutput        = fullModel(cutLine);
            obj.targetFunction   = fullModel;
            obj.nOrder           = nOrder;
            obj.polynomialBasisType = polyType;
            obj.inputSampler = inputSampler;
        end
        
        function obj = run(obj)
            obj = obj.setPolynomials();
            obj = obj.initializeBasisMatrices(obj.nSamples);
            obj = obj.mainEffectSensitivity(1000);
        end
        function obj = initializeBasisMatrices(obj, nSamples)
            %% Project original points onto the polynomial basis for single
            % component HDMR functions
            %% Creates nDim matrices of dimension nSamples x nOrder
            L(obj.nDimensions).matrix = zeros([nSamples, obj.nOrder]);
            L(obj.nDimensions).sampleSet = zeros([nSamples, 1]);
            L(obj.nDimensions).F = zeros([nSamples, 1]);
            cutLineMatrix = repmat(obj.cutLine, nSamples, 1);
            L(obj.nDimensions).C = zeros([obj.nOrder, 1]);
            for iDim = 1:obj.nDimensions
                tempSamples = obj.inputSampler(nSamples);
                L(iDim).sampleSet = tempSamples(:, iDim);
                L(iDim).matrix = zeros([nSamples, obj.nOrder]);
                for iOrder = 0:obj.nOrder - 1
                    index = iOrder + 1;
                    L(iDim).matrix(:, index) = obj.polyFunc(...
                        iOrder, L(iDim).sampleSet);
                end
                tempInputs = cutLineMatrix;
                tempInputs(:, iDim) = L(iDim).sampleSet;
                L(iDim).F = obj.targetFunction(tempInputs) - repmat(obj.cutOutput, nSamples, 1);
                L(iDim).C = obj.hdmr_compressedSensing(...
                    L(iDim).matrix,...
                    L(iDim).F);
            end
            obj.L_struct = L;
        end
        
        function obj = initialize2dBasisMatrices(obj)
           L = obj.L_struct;
            
           [nImportantDim, ~] = size(obj.importantDims);
           %% Compute Psi_2 matrices
           for iDim = 1:nImportantDim
               iSelection = importantDims(iDim);
               for jDim = iDim:nImportantDim
                   jSelection = importantDims(jDim);
                   Psi(iSelection, jSelection).matrix = kron(L(iSelection).matrix, ...
                       L(jSelection).matrix);
                   [mI, ~] = size(L(iSelection).sampleSet);
                   [mJ, ~] = size(L(jSelection).sampleSet);
                   nF = mI*mJ;
                   cutLineMatrix = repmat(obj.cutLine, nSamples, 1);
                   iSamples = reshape(repmat(L(iSelection).sampleSet', mJ, 1), mI*mJ);
                   jSamples = repmat(L(jSelection).sampleSet, mI, 1);
                   cutLineMatrix(:, iSelection) = iSamples;
                   cutLineMatrix(:, jSelection) = jSamples;
                   Psi(iSelection, jSelection).F = obj.targetFunction(cutLineMatrix);
                   Psi(iSelection, jSelection).C = obj.hdmr_compressedSensing(...
                       Psi(iSelection, jSelection).matrix,...
                       Psi(iSelection, jSelection).F);
               end
           end
           
           obj.Psi_struct = Psi;
           
        end
        
        function obj = mainEffectSensitivity(obj, nSamples)
            samples = obj.inputSampler(nSamples);
            disp(size(samples));
            L = obj.L_struct;
            Variances = zeros([1, obj.nDimensions]);
            for iDim = 1:obj.nDimensions
                Psi = zeros(nSamples, obj.nOrder);
                for iOrder = 0:obj.nOrder - 1
                    Psi(:, iOrder + 1) = obj.polyFunc(iOrder, samples(:, iDim));
                end
                f = Psi * L(iDim).C;
                Variances(iDim) = var(f);
            end
            obj.mainEffectSensivitivies = Variances./(sum(Variances));
        end
        
        function C = hdmr_compressedSensing(obj, Psi, f)
            B = f;
            Aineq = [];                                                     % Inequality Constraint matrix
            bineq = [];                                                     % Inequality constraint
            [~, N] = size(Psi);
            lb = zeros([N, 1]);
            ub = [];
            Aeq = Psi;                                                      % Equality Constraint matrix
            beq = B;                                                        % Equlity constraint
            %objective = @(c) sum(abs(c));
            C = linprog(ones([1,N]),Aineq,bineq,Aeq,beq, lb, ub);
        end
        
        function obj = setPolynomials(obj)
            if strcmp(obj.polynomialBasisType, 'bessel')
                obj.polyFunc = @(order, x) besselP(order, x);
            elseif strcmp(obj.polynomialBasisType,'hermite')
                obj.polyFunc = @(order, x) hermiteP(order, x);
            elseif strcmp(obj.polynomialBasisType, 'legendre')
                obj.polyFunc = @(order, x) legendreP(order, x);
            else
                error("Polynomial basis type not supported: see cutHDMR_compressedSensing.m for supported types /n")
            end  
        end
        
        function [sampSorted, fSorted] = plotSingleComponentFunction(obj, iDim, nSamples)
            rng(0);
            n = get(gcf,'Number');
            figure(n + 1)
            samples = obj.inputSampler(nSamples);
            disp(size(samples));
            L = obj.L_struct;

            Psi = zeros(nSamples, obj.nOrder);
            for iOrder = 0:obj.nOrder - 1
                Psi(:, iOrder + 1) = obj.polyFunc(iOrder, samples(:, iDim));
            end
            f = Psi * L(iDim).C + repmat(obj.cutOutput, nSamples, 1);
            [sampSorted, fSorted] = sortTwoArrays(samples(:, iDim), f);
            plot(sampSorted, fSorted);
        end
    end
    
end