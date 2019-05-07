% Earlier version of this function (version 1).
% Replaced on 8/14/2007.
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
function [c,A] = AffineEmaxBreakpointsV1(a,b)
    % Preallocate c for speed.
    M = length(a);
    c = zeros(1,M+1);
    
    % Step 0
    i = 1;
    c(1) = -inf;
    c(2) = +inf;
    A = [1];
    
    while(i<M)
        % Step 1
        k = i+1;
        c(1+i+1) = +inf;
        % c(1) = -inf; % Always true.
        A = [A, k];
        kindex = length(A);
    
        % Step 2
        % A has at least two entries.  k is the last entry, and j is the
        % second to last.
        j = A(kindex-1);
    
        while (j>0)
            changed_j = 0;
            while (1)
                % At the beginning of this loop, we will always have j>0
                % since, on the previous loop, we would have exited if l=0
                % since c(1+l)=-inf, and c(1+j) is a finite number and
                % hence strictly larger than -inf.
                
                % Step 3
                if (kindex == 2), l=0; else, l=A(kindex-2); end
    
                % Step 4
                c(1+j) = (a(j) - a(k))/(b(k)-b(j));
    
                % Step 5
                if (c(1+j) > c(1+l)), break; end
                % Discard j from A.  j is the element in A before k.
                A = A([[1:kindex-2],[kindex:length(A)]]);
                %A = setdiff(A,[j]);
                % Since A has one fewer elements, decrement kindex.
                kindex = kindex - 1; 
                % After j's removal, l is the second before k in A,
                % so set j to l to restore j.  l will be reset to the
                % element in A before j, and second before k.
                j=l; % equiv to j=A(kindex-1);
                changed_j = 1;
            end
    
            if (changed_j == 0), break; end
            k=j; kindex = kindex-1;
            j=l;
        end
    
        i=i+1;
    end
    
end