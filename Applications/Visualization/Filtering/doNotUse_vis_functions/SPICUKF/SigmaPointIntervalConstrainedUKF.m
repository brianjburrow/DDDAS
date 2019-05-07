function input = SigmaPointIntervalConstrainedUKF(input)
    %% SigmaPointIntervalConstrainedUKF(input)
    % Input
    % input.x = state vector (column vector N x 1)
    x = input.x;
    P = input.P;
    lambda = input.lambda;
    n = input.n;
    F = input.F;
    G = input.G;
    data = input.data;
    Q = input.Q;
    u = input.u;
    R = input.R;
    lb = input.lb;
    ub = input.ub;
    x0 = x;
    try
        % if ukf == 1, then perform UKF w/o bounds
        ukf = input.ukf;
    catch
        ukf = 0;
    end

    %% Unpack input object
    %% project state estimate onto bounds
    
    %% Forecast
    [X, weights] = SPICUKF_sigma_ICUT(x, P, lb, ub, n, lambda);
    % Pass sigma points through data evolution model  
    for idx = 1:(2*n+1)
        X(:,idx) = F(X(:,idx), u);
    end
    % Forcast updated mean
    x = zeros(length(x), 1);
    for idx = 1:(2*n+1)
        x = x + X(:,idx)*weights(idx);
    end
    

    % Forecast updated covariance
    P = zeros(size(P));
    for idx = 1:(2*n)
        P = P + weights(idx) * (...
            (X(:,idx) - x) *...
            (X(:,idx) - x)'...
            );
    end

    P = P + Q;
    % update weights and sigma points
    [X, weights] = SPICUKF_sigma_ICUT(x, P, lb, ub, n, lambda);
    
    %% Data Assimilation
    if ~ ukf
        % precalc R^-1
        R_inv = R\eye(size(R));
        % precalc P^-1
        P_inv = P\eye(size(P));
        X_prev = X;
        X = zeros(size(X));

        % optimization options
        try
            A = input.A;
            b = input.b;
        catch
            A = [];
            b = [];
        end
        Aeq = [];
        beq = [];
        nonlcon = [];
        opts1 = optimoptions('fmincon');
        %opts1.PlotFcn = 'optimplotstepsize';
        %opts1.Algorithm = 'sqp-legacy';
        opts1.Display = 'off';
        dataFunc_solo = @(xx) G(xx, u);
        for idx = 1:2*n + 1
            Xtemp = zeros(size(X));
            fval = zeros(1, 2*n+1);
            for dmx = 1:2*n + 1
                % typically use X_prev(:,idx) as initial point.  Attempting x
                clear func
                func = @(xx) J_2(xx, X_prev(:,idx),...
                    data, dataFunc_solo, R_inv, P_inv);
                [Xtemp(:, dmx), fval(dmx)] = fmincon(...
                    func, X_prev(:,dmx), A, b, Aeq, beq, lb, ub, nonlcon, opts1);
            end
            [~, minIdx] = min(fval);
            X(:,idx) = Xtemp(:, minIdx);
        end
        disp("Assimilating Data: new sigma points")
        disp(X)
    else
        % Predict measurements
        Y = zeros(length(data), 2*n+1);
        for idx = 1:(2*n+1)
            Y(:, idx) = G(X(:,idx), u);
        end
        y = zeros(length(data), 1);
        for idx = 1:(2*n+1)
            y = y + Y(:,idx)*weights(idx);
        end
        P_yy = zeros(size(R));
        P_xy = zeros(length(P), length(R));

        for idx = 1:(2*n)
            P_yy = P_yy + weights(idx) * (...
                (Y(:,idx) - y) *...
                (Y(:,idx) - y)'...
                );

            P_xy = P_xy + weights(idx) * (...
                (X(:,idx) - x) *...
                (Y(:,idx) - y)'...
                );
        end
        P_yy = P_yy + R;
        %% Calculate Kalman Gain
        K = P_xy * (P_yy\eye(size(P_yy)));
        x = x + K*(data - y);
        P = P - K*P_yy*K';
    end
    % Update mean and covariance

    x = zeros(length(x), 1);
    for idx = 1:(2*n+1)
        x = x + X(:,idx)*weights(idx);
    end
    P = zeros(size(P));

    for idx = 1:(2*n)
        P = P + weights(idx)' * (...
            (X(:,idx) - x) *...
            (X(:,idx) - x)'...
            );
    end
    %% Repackage input

    input.x = x;
    % add small diagonal constant to avoid ill conditioning, b/c P can
    % become non-positive semi-definite, causing Chol to fail
    input.P = P + 1e-14 * eye(size(P));
    %input.P = P + Q;
    input.X = X;
    function objective = J_2(X_column, X_prev_column, data, ...
            dataFunc,R_inv, P_inv)
        %minim = min(data);
        %maxim = max(data);
        %data = (data - minim)./(maxim - minim);
        %scaledDataFunc = (dataFunc(X_column) - minim)./(maxim - minim);
        scaledDataFunc = dataFunc(X_column);

        objective = (data - scaledDataFunc)'*(R_inv)*...
            (data - scaledDataFunc) +...
            (X_column - X_prev_column)' *...
            (P_inv)*...
            (X_column - X_prev_column);
    end
    
%     function objective = J_2(X_column, X_prev_column, data, ...
%             dataFunc,R_inv, P_inv, minim, maxim)
%         scaledData = (dataFunc(X_column) - minim)./(maxim - minim);
%         objective = (data - dataFunc(X_column))'*(R_inv)*...
%             (data - dataFunc(X_column)) +...
%             (X_column - X_prev_column)' *...
%             (P_inv)*...
%             (X_column - X_prev_column);
%     end
end