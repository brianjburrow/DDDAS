% logy = LogEmaxAffine(a,b)
% Calculates log(Exp[max_x a_x + b_x Z]-max_x a_x), where Z is a standard
% normal random variable and a,b are 1xM input vectors.
function [logy, a,b,c] = LogEmaxAffine(a,b)
    if (any(isnan(a)) || any(isnan(b)))
        warning('a or b is NaN');
    end
    assert(all(isreal(a)));
    assert(all(isreal(b)));

    a = ConvertToColumnVector(a);
    b = ConvertToColumnVector(b);

    % Check that a and b are column vectors of the right size
    sa = size(a);
    if (sa(2) ~= 1 || any(size(a) ~= size(b)))
        error('LogEmaxAffine: a and b must be column vectors of the same size');
    end
    
    [a,b] = AffineBreakpointsPrep(a,b);
    if (length(a)==1)
        logy = -Inf;
        return
    end
    [c, keep] = AffineBreakpoints(a,b); 
    
    a = a(keep);
    b = b(keep);
    c = c([1,keep+1]);
    M = length(keep);
    assert(all(isreal(c)));

    % I need logbdiff=log(diff(b)).  I thought that the following code would be
    % more numerically stable, able for example to distinguish cases like 
    % logb = [-25 -.3] vs. logb = [-35 -.3], but it doesn't seem to be able to.
    % Indeed, in the debugging output that I have below, the difference was 0.
    %{
    logb = log(abs(b)); % If b is 0, this is -Inf.
    sgnb = sign(b); % If b is 0, this is 0.
    logbdiff = zeros(size(c(2:M)));
    for i=1:length(b)-1
	[logbdiff(i),logbdiffsgn] = LogPlusExpSigned(logb(i),sgnb(i),logb(i+1),-sgnb(i+1));
	%assert(logbdiffsgn>=0);  % The b are distinct, so bdiff(i) can't be 0.
    end
    disp(sprintf('log(b)=%s log(diff(b))=%g logbdiff=%g difference=%g',mat2str(log(b)),log(diff(b)),logbdiff,log(diff(b))-logbdiff));
    %}
    logbdiff = log(diff(b))';   
    assert(M>=2);
    logy = LogSumExp(logbdiff+LogEI(-abs(c(2:M))));
    assert(isreal(logy));
end
