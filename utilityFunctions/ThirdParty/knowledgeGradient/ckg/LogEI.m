% Returns the log of Exp[(s+Z)^+], where s is a constant and Z is a standard
% normal random variable.  For large negative arguments Exp[(s+Z)^+] function
% is close to 0.  For large positive arguments, the function is close to the
% argument.  For s large enough, s>-10, we use the formula
% Exp[(s+Z)^+] = s*normcdf(s) + normpdf(s).  For smaller s we use an asymptotic
% approximation based on Mill's ratio.  EI stands for "expected improvement",
% since Exp[(s+Z)^+] would be the log of the expected improvement by measuring
% an alternative with excess predictive mean s over the best other measured
% alternative, and predictive variance 0.
function logy = LogEI(s)

% Use the asymptotic approximation for these large negative s.  The
% approximation is derived via:
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-|s|normcdf(-|s|)/normpdf(s)]
% and noting that normcdf(-|s|)/normpdf(s) is the Mill's ratio at |s|, which is
% asymptotically approximated by |s|/(s^2+1) [Gordon 1941, also documented in
% Frazier,Powell,Dayanik 2009 on page 14].  This gives,
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-s^2/(s^2+1)] = normpdf(s)/(s^2+1).
i=find(s<-10);
if (length(i)>0)
    logy(i) = LogNormPDF(s(i)) - log(s(i).^2 + 1);
end

% Use straightforward routines for s in the more numerically stable region.
i=find(s>=-10);
if (length(i)>0)
    logy(i) = log(s(i).*normcdf(s(i))+normpdf(s(i)));
end

assert(all(isreal(logy)));
