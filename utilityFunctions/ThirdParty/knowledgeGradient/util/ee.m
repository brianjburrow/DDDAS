% Creates an nxn matrix with a single 1 at entry mm, and 0s elsewhere.
function y = ee(m,n)
    y = zeros(n);
    y(m,m) = 1;
end