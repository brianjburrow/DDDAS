function theta = SPICUKF_fill_theta(n, lambda, e, d, x, S)
    %% Inputs
    % n : number of dimensions (scalar)
    % lambda: scaling parameter (scalar)
    % e: upper bound (n x 1) array
    % d: lower bound (n x 1) array
    % x: current state (n x 1) array
    % S: [P -P] where P is the state covariance
    % S: if P is (n x n) array, then S is (n x 2n) array
    %% outputs
    % THETA vector necessary for SPICUKF algorithm
    % See 'Robust and reliable estimation via Unscented Recursive
    % Nonlinear Dynamic Data Reconciliation' by Vachhani
    % Be wary of 'On unscented Kalman filtering with state interval
    % constraints' by Teixeira.  I'm not sure if different, or I can't
    % code...
    baseVal = sqrt(n + lambda);
    S_siz = size(S);
    thetaM = zeros(S_siz);
    
    for idx = 1:S_siz(1)
        for dmx = 1:S_siz(2)
            if S(idx, dmx) == 0
                thetaM(idx, dmx) = baseVal;
            elseif S(idx, dmx) > 0
                thetaM(idx, dmx) = min(baseVal, ...
                    (e(idx) - x(idx))/S(idx, dmx));
            else
                thetaM(idx, dmx) = min(baseVal, ...
                    (d(idx) - x(idx))/S(idx, dmx));
            end
        end
    end

    theta = zeros(2*n, 1);
    for idx = 1:2*n
        theta(idx) = min(thetaM(:,idx));
    end
end