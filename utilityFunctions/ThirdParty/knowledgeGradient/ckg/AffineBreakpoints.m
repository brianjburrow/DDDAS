% Inputs are two M-vectors, a and b.
% Requires that the b vector is sorted in increasing order.
% Also requires that the elements of b all be unique.
% This function is used in AffineEmax, and the preparation of generic
% vectors a and b to satisfy the input requirements of this function are
% shown there.
%
% The output is an (M+1)-vector c and a vector A ("A" is for accept).  Think of
% A as a set which is a subset of {1,...,M}.  This output has the property
% that, for any i in {1,...,M} and any real number z,
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
function [c,A] = AffineBreakpoints(a,b)
    % Preallocate for speed.  Instead of resizing the array A whenever we add
    % to it or delete from it, we keep it the maximal size, and keep a length
    % indicator Alen telling us how many of its entries are good.  When the
    % function ends, we remove the unused elements from A before passing
    % it.
    M = length(a);
    c = zeros(1,M+1);
    A = zeros(1,M);
    
    % Step 0
    i=0;
    c(1+i) = -inf;
    c(1+i+1) = +inf;
    A(1) = 1;
    Alen = 1;
    
    for i=[1:M-1]
        c(1+i+1) = +inf;
        while(1)
            j = A(Alen); % jindex = Alen
            c(1+j) = (a(j) - a(i+1))/(b(i+1)-b(j));
	    % The if statement below replaces these lines from version 2 of the
	    % function.
	    %    kindex = jindex-1 = Alen-1
            %    if kindex > 0 && c(1+j)<=c(1+A(kindex))
            if Alen > 1 && c(1+j)<=c(1+A(Alen-1))
		Alen = Alen-1; % Remove last element j
                % continue in while(1) loop
            else
                break % quit while(1) loop
            end
        end
	A(Alen+1) = i+1;
	Alen = Alen + 1;
    end
    A = A(1:Alen);
