% Tests NoninformativeBayesUpdate and NoninformativeBayesUpdateSqrt on a
% problem where the alternatives are numbers between 1 and 15, and an
% alternative has features corresponding to which bits are set in its number.
% So alternative 7 has features 1,2, and 3 since 7 = 1110 in binary.
% This function goes through 4 measurements and uses these two update functions
% to get the posterior.  It compares the final posterior means and marginal
% variances to what they should be theoretically and reports any differences.
% These differences should be small, on the order of machine precision.  In my
% tests on Nov 12, 2008 I see they are less than 6.3e-14.
function TestNoninformativeBayesUpdate()
for i=1:15
for j=1:15
	bits_in_common = bitand(i,j);
	num_in_common = sum(bitget(bits_in_common,1:8));
	B(i,j) = num_in_common;
end
end
% This covariance matrix I have just constructed is not quite positive definite
% due to matlab rounding errors.  Fix this with a fudge factor:
%B = B-min(eig(B))*eye(15);
%sqrtB = chol(B,'lower')
sqrtB = cholinc(sparse(B),1e-15)';
 
A = zeros(15);
A_withsqrt = A;
mu = zeros(15,1);
mu_withsqrt = mu;

x = [1,2,4,8];
y = [0.8147,   0.9058,    0.1270,    0.9134];
lambda = 1;

for n=1:4
	[mu,A,B] = NoninformativeBayesUpdate(mu,A,B,x(n),lambda,y(n));
	[mu_withsqrt,A_withsqrt,sqrtB] = NoninformativeBayesUpdateSqrt(mu_withsqrt,A_withsqrt,sqrtB,x(n),lambda,y(n));
	assert(all(all(sqrtB==tril(sqrtB))));
end

% Calculate theoretically what the posterior ought to be.  We measured each
% substituent once, and they are independent of each other with variance 1, so
% the posterior on each alternative should be equal to the sum of the y for the
% bits it has set, and its variance should be equal to the number of bits it
% has set.
for n=1:15
	mu_prediction(n) = 0;
	var_prediction(n) = 0;
	for bit=1:4
		if (bitget(n,bit))
			mu_prediction(n) = mu_prediction(n) + y(bit);
			var_prediction(n) = var_prediction(n) + 1;
		end
	end
end


disp('Column 1=NoninformativeBayesUpdate');
disp('Column 2=NoninformativeBayesUpdateSqrt');
disp('Column 3=Theoretical');
disp('posterior means');
disp([mu,mu_withsqrt,mu_prediction']);
disp('posterior marginal variances')
disp([diag(A),diag(A_withsqrt),var_prediction']);

disp(sprintf('discrepancy of %g between NoninformativeBayesUpdateSqrt and the theoretical posterior mean', max(abs(mu_withsqrt'-mu_prediction))));
disp(sprintf('discrepancy of %g between NoninformativeBayesUpdateSqrt and NoninformativeBayesUpdate posterior mean', max(abs(mu_withsqrt'-mu'))));

disp(sprintf('discrepancy of %g between NoninformativeBayesUpdateSqrt and the theoretical posterior variance', max(abs(diag(A_withsqrt)-var_prediction'))));
disp(sprintf('discrepancy of %g between NoninformativeBayesUpdateSqrt and NoninformativeBayesUpdate posterior variance', max(abs(diag(A_withsqrt)-diag(A)))));
