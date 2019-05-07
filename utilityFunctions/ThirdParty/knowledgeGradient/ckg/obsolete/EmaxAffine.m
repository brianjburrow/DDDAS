% Calculates \Exp\max_x a_x + b_x Z, where Z is a standard normal random
% variable and a,b are 1xM input vectors.  It uses the algorithm and some
% of the notation from the Correlated Normal Knowledge-gradient paper,
% pages 6-7 in the 07062007 draft.
function emax = EmaxAffine(a,b)
    if (any(isnan(a)) || any(isnan(b)))
	%warning('a or b is NaN');
    end
    assert(all(isreal(a)));
    assert(all(isreal(b)));

    % Check that a and b are column vectors of the right size
    sa = size(a);
    if (sa(2) ~= 1 || any(size(a) ~= size(b)))
        error('AffineEmax: a and b must be column vectors of the same size');
    end
    
    [a,b] = AffineBreakpointsPrep(a,b);
    [c, keep] = AffineBreakpoints(a,b);
    
    a = a(keep);
    b = b(keep);
    c = c([1,keep+1]);
    M = length(keep);
    
    assert(all(isreal(c)));
    term1 =  normcdf(c(2:M+1)) - normcdf(c(1:M));
    term2 = -normpdf(c(2:M+1)) + normpdf(c(1:M));
    emax = term1*a + term2*b;

    if (emax == max(a) && any(b~=0))
	disp(sprintf('We computed Exp[max_x a_x + b_x Z] = max_x a_x even though b is not identically 0.  a=%s b=%s c=%s',mat2str(a),mat2str(b),mat2str(c)));
	
    end
end
