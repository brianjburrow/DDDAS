function z = ConvertRotatedDomain(x, lowBound, upBound)
    % Convert points in an arbitrary D-dimensional hypercube to the [-1, 1]^D
    % x is an nSamples x D matrix;

    [nSamples, D] = size(x);
    lowBound = reshape(lowBound, [1, D]);
    lowBound = repmat(lowBound, [nSamples, 1]);

    upBound = reshape(upBound, [1, D]);
    upBound = repmat(upBound, [nSamples, 1]);

    z = (2*x - (upBound + lowBound))./(upBound - lowBound);
end
