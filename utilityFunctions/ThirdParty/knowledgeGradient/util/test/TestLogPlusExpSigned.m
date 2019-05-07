function ok = TestLogPlusExpSigned()
	ok = 1;
	difference  = zeros(1,100);
	for n=1:100
		[difference(n),this_ok]=helper();
		ok = this_ok && ok;
	end
	%hist(difference); % can be useful for debugging.
	if (ok)
		disp('TestLogPlusExpSigned: OK');
	else
		disp('TestLogPlusExpSigned: FAILED');
	end
end

function [difference, ok] = helper()
	a = normrnd(0,1);
	b = normrnd(0,1);
	c = a+b;
	%logc = log(abs(c));
	%sgnc = sign(c);
	[logd sgnd]=LogPlusExpSigned(log(abs(a)),sign(a),log(abs(b)),sign(b));
	d = sgnd*exp(logd);
	difference = c-d;
	if (abs(difference)>1e-8)
		ok = 0;
		disp(sprintf('Error,a=%g b=%g logd=%g sgnd=%g',a,b,logd,sgnd));
	end
	ok = 1;
end
