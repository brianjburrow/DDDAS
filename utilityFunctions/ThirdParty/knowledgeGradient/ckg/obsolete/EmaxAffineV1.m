% Calculates \Exp\max_x a_x + b_x Z, where Z is a standard normal random
% variable and a,b are 1xM input vectors.  It uses the algorithm and some
% of the notation from the Correlated Normal Knowledge-gradient paper,
% pages 6-7 in the 07062007 draft.
function emax = AffineEmax(a,b)
    M = length(a);
    % Form a matrix for which ba(x,1) is the slope b(x) and ba(x,2) is the
    % y-intercept a(x).  Sort this matrix in ascending order of slope,
    % breaking ties with the y-intercept.
    ba = [b, a];
    ba = sortrows(ba,[1,2]);
    a = ba(:,2);
    b = ba(:,1);
    
    % Need to discard those entries with the same slope!!
    
    % Algorithm for computing c.
    
    % Step 0
    i = 1;
    oldc= [-inf, +inf];
    A = [1];
    
    while(i<M)
        % Step 1
        k = i+1;
        c(1+i+1) = +inf;
        c(1) = -inf;
        A = [A, k];
        kindex = length(A);
    
        % Step 2
        jindex = kindex-1;
        j = A(jindex);
    
        while (j>0)
            while (1)
                % Step 3
                lindex = jindex - 1;
                if (lindex == 0), l=0; else, l=A(lindex); end
    
                % Step 4
                c(1+j) = (a(j) - a(k))/(b(k)-b(j));
    
                % Step 5
                if (c(1+j) > oldc(1+l)), break; end
                A = A([[1:jindex-1],kindex]); % discard j.
                kindex = kindex-1;
                j=l; jindex = lindex;
                % Go to step 3.
            end
    
            k=j; kindex = jindex;
            j=l; jindex = lindex;
        end
    
        i=i+1;
        oldc = c;
    end
    
    
    a = a(A);
    b = b(A);
    M = length(A);   
    
    term1 =  normcdf(c(2:M+1)) - normcdf(c(1:M));
    term2 = -normpdf(c(2:M+1)) + normpdf(c(1:M));
    emax = sum(a' .* term1 + b' .* term2);
end