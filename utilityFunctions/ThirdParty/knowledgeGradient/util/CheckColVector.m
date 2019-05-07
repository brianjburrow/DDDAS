% Checks if a vector is a column vector, and if it isn't, produces an error.
% This is for use in some code in place of ConvertToColumnVector, if converting
% the vector without telling the user would be dangerous.
function CheckColVector(x)
[nrows,ncols] = size(x);
if (nrows == 1 && ncols > 1)
	error('Passed vector is not a column vector');
end
