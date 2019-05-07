classdef LegendrePolynomialCurve < handle
    % Fits an n^th order, univariate legendre polynomial through a set of 
    % training data.  
    properties
        order
        designMatrix
        outputs
        nSamples
        nOutputDimensions
        nInputDimensions
        coefficientMatrix    %
        regressionCoefficients
    end
    
    methods
        function obj = LegendrePolynomialCurve(order, design, outputs)
            %% Inputs
            % order: Polynomial order of the curve fitting                  Scalar
            % design: Design matrix (i.e., input training data)             nSamples x 1
            % outputs: Outputs corresponding to the input training data     nSamples x nOutputDimension  (i.e., fit nOutputDim independent polynomials)
            obj.order = order;
            obj.designMatrix = design;
            obj.outputs = outputs;
            [obj.nSamples, obj.nOutputDimensions] = size(design);
            [~, obj.nInputDimensions] = size(outputs);
            obj.fit()
        end
        
        function obj = fit(obj)
            disp('Fitting polynomial coefficients')
            obj.coefficientMatrix = zeros(obj.nSamples, obj.order + 1);
            for iOrder = 1:obj.order + 1
                obj.coefficientMatrix(:, iOrder) = obj.evalPolynomials(...
                    obj.designMatrix, iOrder - 1 ...
                );                                                           % Must subtract 1 so that order starts with 0, but matlab indexs from 1, ...
            end
            obj.regressionCoefficients = mldivide(...
                obj.coefficientMatrix, obj.outputs...
                );
        end
        function predictions = predict(obj, X)
            [nSamp, ~] = size(X); 
            coefMat = zeros(nSamp, obj.order + 1);
            for iOrder = 1:obj.order + 1
                coefMat(:, iOrder) = obj.evalPolynomials(X, iOrder - 1);
            end
            predictions = coefMat * obj.regressionCoefficients;
        end
        
        function output = evalPolynomials(obj, X, order)
            if order == 0
                output = ones(size(X));
            elseif order == 1
                output = X;
            elseif order == 2
                output = 0.5 * (3*X.^2 - 1);
            elseif order == 3
                output = 0.5 * (5*X.^3 - 3*X);
            elseif order == 4
                output = 0.125 * (35 * X.^4 - 30.*X.^2 + 3);
            elseif order == 5
                output = 0.125 * (63 * X.^5 - 70 * X.^3 + 15 * X);
            else
                error("Legendre Polynomials of 6th order or higher not implemented")
            end
        end
    end
end