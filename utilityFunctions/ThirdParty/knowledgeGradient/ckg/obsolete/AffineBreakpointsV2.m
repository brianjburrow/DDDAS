% Earlier version of this function (version 2).
% Replaced on 5/30/2008.
%
% Inputs are two M-vectors, a and b.
% Requires that the b vector is sorted in increasing order.
% Also requires that the elements of b all be unique.
% This function is used in AffineEmax, and the preparation of generic
% vectors a and b to satisfy the input requirements of this function are
% shown there.
%
% The output is an (M+1)-vector c and a vector A ("A" is for accept).
% Think of A as actually a set, actually a subset of {1,...,M}.
% This output has the property that, for any i in {1,...,M} and any real
% number z,
%   i \in argmax_j a_j + b_j z
% iff
%   i \in A and z \in [c(j+1),c(i+1)],
%   where j = sup {0,1,...,i-1} \cap A.
%
% A note about indexing:
% Since Matlab does not allow indexing from 0, but instead requires
% indexing from 1, what is called c_i in the paper is written in matlab as
% c(1+i).  This is because in the paper we reference c_0.  For the vectors
% a and b, however, we don't need to reference a_0 or b_0, so we reference
% a_i and b_i by a(i) and b(i) respectively, rather than a(i+1) or b(i+1).
% 
function [c,A] = AffineEmaxBreakpoints(a,b)
    % Preallocate c for speed.
    M = length(a);
    c = zeros(1,M+1);
    
    % Step 0
    i=0;
    c(1+i) = -inf;
    c(1+i+1) = +inf;
    A = [1];
    
    for i=[1:M-1]
        c(1+i+1) = +inf;
        while(1)
            jindex = length(A);
            j = A(jindex);
            c(1+j) = (a(j) - a(i+1))/(b(i+1)-b(j));
            kindex = jindex-1;
            if kindex > 0 && c(1+j)<=c(1+A(kindex))
                A = A(1:kindex); % Remove last element j
                % continue in while(1) loop
            else
                break % quit while(1) loop
            end
        end
        A = [A,i+1];
    end
