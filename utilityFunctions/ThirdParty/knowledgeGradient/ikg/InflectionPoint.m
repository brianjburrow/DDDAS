% Finds the inflection point of the function that returns the KG factor as a
% function of the number of measurements, n.  The inflection point is the point
% where the function has 0 second derivative, and switches from being convex to
% being concave.
function n = InflectionPoint(delta, beliefvar, noisevar)
term1 = delta^4 + 14*beliefvar*delta^2+beliefvar^2;
term2 = delta^2 - beliefvar + sqrt(term1);
n = term2 * noisevar / (8*beliefvar^2);
