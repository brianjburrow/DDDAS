function [newLowBounds, newUpBounds] = FindRotatedDomain(...
    lowerBounds, ...
    upperBounds,...
    rotationMatrix)
%% FindRotatedDomain
% Finds min and maximum of an input domain after the basis is rotated.
% Takes in a vector of lower bounds and a vector of upper bounds, where each
% is an nDim x 1 vector.
% Rotation matrix must be a rotation matrix of nDim x nDim dimensions.
%% How it works
% This function takes the lower and upperbounds and finds all vertices of
% the hypercube defined by the bounds.  Then it rotates each of those vertices
% and computes the min and max ranges after that rotation.  This does not
% mean that the function is defined over this entire domain.

nDim = length(lowerBounds);

if nDim ~= length(upperBounds)
    error("Dimensions of the lower bound and upper bound vectors must be the same")
end

bounds = [lowerBounds, upperBounds];
design = fullfact(2*ones(1, nDim))';

nDesign = length(design(1, :));
vertices = zeros(nDim, nDesign);
for iDesign = 1:nDesign
    for iDim = 1:nDim
        vertices(iDim, iDesign) = bounds(iDim, design(iDim, iDesign));
    end
end

disp(vertices)

vertices = rotationMatrix * vertices;

newLowBounds = min(vertices, [], 2);
newUpBounds = max(vertices, [], 2);
