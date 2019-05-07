% Creates a nx1 matrix with a single 1 at entry m, and 0s elsewhere.
function y = e(m,n)
    y = zeros(n,1);
    y(m) = 1;
end