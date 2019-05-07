% Let Y= log(sgna*exp(a)+sgnb*exp(b)), where sgn[ab] are -1,0,+1.  Then
% y=log(abs(Y)) and sgny=sign(Y).  This is computed in a numerically stable
% way.
function [y,sgny]=LogPlusExpSigned(a,sgna,b,sgnb)
	if (sgna==0)
		y = b;
		sgny = sgnb;
	elseif (sgnb==0)
		y = a;
		sgny = sgna;
	elseif (sgna == sgnb)
		sgny = sgna;
		y = LogSumExp([a b]);
	elseif (a==b) % and signs are opposite
		y = -inf;
		sgny = 0;
	elseif (a>b)
		sgny = sgna;
		y = LogMinusExp(a,b);
	else
		sgny = sgnb;
		y = LogMinusExp(b,a);
	end
end
