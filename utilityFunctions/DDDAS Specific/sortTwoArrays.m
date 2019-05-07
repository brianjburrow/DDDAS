function [a_sorted, newB] = sortTwoArrays(A, B)
    [a_sorted, a_order] = sort(A);
    newB = B(a_order,:);
end