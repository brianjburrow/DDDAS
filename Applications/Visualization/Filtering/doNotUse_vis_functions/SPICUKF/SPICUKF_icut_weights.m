function weights = SPICUKF_icut_weights(lambda, sum_theta_vec, n, theta_vec)

        a = (2*lambda - 1) / ...
            (...
            2*(n + lambda)*...
            (sum_theta_vec - (2*n+1)*sqrt(n+lambda))...
            );
        
        b = (1/(2*(n + lambda))) - ...
            (2*lambda - 1) / ...
            (...
            (2*sqrt(n + lambda)) * (...
                sum_theta_vec - (2*n + 1)*sqrt(n + lambda)...
                )...
            );
        
        weights = zeros(2*n+1,1);
        weights(1) = b;
        weights(2:end) = a*theta_vec + b;
    end