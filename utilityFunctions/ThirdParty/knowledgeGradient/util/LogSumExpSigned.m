% Let Y= log(\sum_i sgna_i*exp(a_i)), where sgn[a_i] is -1,0,+1.  Then
% y=log(abs(Y)) and sgny=sign(Y).  This is computed in a numerically stable
% way.
function [y,sgny]=LogSumExpSigned(a,sgna)
assert(length(a) == length(sgna));
y = a(1);
sgny = sgna(1);
for i=2:length(a)
	[y,sgny] = LogPlusExpSigned(y,sgny,a(i),sgna(i));
end
