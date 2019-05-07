function [X, weights] = SPICUKF_sigma_ICUT(x, P, d, e, n, lambda)
    %% Inputs
    % x : current state estimate (n x 1) array
    % P : current state covariange (n x n) array
    % d : lower bound on state values (n x 1) array
    % e : upper bound on state values (n x 1) array
    % n : number of dimensions (scalar)
    % lambda: scaling parameter (scalar)
    %% Returns
    % X : set of sigma points for Sigma Point Interval Constrained UKF
    % weights: corresponding weights for each sigma point

        a = chol(P);
        S = [a, -a];
        P_dubs = [a, a];
        % find the min of each column of THETA
        thetas = SPICUKF_fill_theta(n, lambda, e, d, x, S);
        sum_thetas = sum(thetas);
        weights = SPICUKF_icut_weights(lambda, sum_thetas, n, thetas);
        
        thetas(n + 1:end) = -1*thetas(n + 1:end);
        X = repmat(x,1, (2*n+1));
        X_old = X;
        aa = zeros(3, length(thetas));
        for idx = 1:2*n
            aa(:,idx) = thetas(idx)*P_dubs(:,idx);
            X(:,idx+1) = X(:,idx+1) + thetas(idx)*P_dubs(:,idx);
        end
    end