function gpObj = appendToMatlabGP(gpObj, X, Y)
    %% NOTE: NOT FUNCTIONAL, I don't know how to update the alpha vector just yet
    %% This function takes in an already trained gp object, and modifies its
    % properties to allow adding new training data dynamically.
    [nAddSamples, ~] = size(X);
    nOrigSamples = gpObj.NumObservations;
    gpObj.NumObservations = nOrigSamples + nAddSamples;
    gpObj.isActiveSetVector = [gpObj.isActiveSetVector; ones([nAddSamples, 1])];   % Update the active set list
    gpObj.Y = [gpObj.Y; Y];
    gpObj.X = [gpObj.X; X];
    gpObj.W = [gpObj.W; gpObj.W(1, 1)*ones([nAddSamples, 1])];
    gpObj.ActiveSetVectors = gpObj.X;                                       % Assumes you trained on the full set
end