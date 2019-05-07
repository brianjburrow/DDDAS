% y=LogMinusExp(a,b)
% Computes y=log(exp(a)-exp(b)) in a numerically stable manner.  Requires that a>=b.
function y=LogMinusExp(a,b)
	assert(a>=b);
	if (a==b)
		y = -Inf;
		return;
	end
	y = a+log1p(-exp(b-a));
end
