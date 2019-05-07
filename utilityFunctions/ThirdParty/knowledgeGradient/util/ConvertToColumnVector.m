% As its name indicates, converts the passed vector x (which may be a row
% vector) to a column vector.  It also produces an error message if a matrix is
% passed, so it may be used to do error checking of inputs.
function y=ConvertToColumnVector(x)
[rows,cols] = size(x);
if (rows > 1 && cols > 1)
  error('Must pass a vector, not a matrix.');
end
if (rows>cols)
  y = x;
else
  y = x';
end
