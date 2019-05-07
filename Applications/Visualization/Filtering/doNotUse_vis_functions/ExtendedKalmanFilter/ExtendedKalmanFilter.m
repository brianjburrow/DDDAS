function [x, P] = ExtendedKalmanFilter(x, P, data, f, h, u, ...
    approxDerivs, F, H, approxStep, Q, R)
%% Performs Extended Kalman Filter (first order)
% Inputs:
%   x is the previous state estimate
%   P is the previous covariance matrix estimate
%   data is data collected at the current time step
%   f is the state transition function.  Must be passed as a function
%   u is the control variables.
% Outputs:
%   x is the updated estimate of the state
%   P is the updated estimate of the covariance matrix
%   h is the data function

%% Calculate derivatives
if approxDerivs
    % approximate derivative using central difference
    % F 
    F = centralDifferenceF(f, [x, u*ones(size(x))], approxStep); 
    % H
    H = centralDifferenceH(h, [x, u*ones(size(x))], approxStep);
end

%% Predict Step
x = f(x, u);

P = F*P*F' + Q;

%% Update Step
%% compute innovation
y = data - h(x, u);

% Calculate S
S = H*P*H' + R;
% Calculate Kalman Gain
K = P*H'*inv(S);

% Update the mean
x = x + K*y;
%Update the variance
P = (eye(size(K*H)) - K*H)*P;
end

function slope = centralDifferenceF(func, centerPoint, step)
    numDim = 4;
    slope = zeros(4, 4);
    for idx = 1:numDim
        vector = zeros(size(centerPoint(:,1)));
        vector(idx) = step;
        slope(idx, :) = (func(centerPoint(:, 1) + vector, centerPoint(:, 2)) - ...
            func(centerPoint(:, 1) - vector, centerPoint(:, 2)))';
        slope(idx, :) = slope(idx, :)./(2*step);
    end
    slope = slope';
end

function slope = centralDifferenceH(func, centerPoint, step)
    numDim = 4;
    slope = zeros(4, 13);
    for idx = 1:numDim
        vector = zeros(size(centerPoint(:,1)));
        vector(idx) = step;
        slope(idx, :) = (func(centerPoint(:, 1) + vector, centerPoint(:, 2)) - ...
            func(centerPoint(:, 1) - vector, centerPoint(:, 2)))';
        slope(idx, :) = slope(idx, :)./(2*step);
    end
    slope = slope';
end