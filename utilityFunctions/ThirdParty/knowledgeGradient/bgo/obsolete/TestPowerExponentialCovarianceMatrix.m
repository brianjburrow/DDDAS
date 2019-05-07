function TestPowerExponentialCovarianceMatrix()
	test1();
	test2();
end
	

function test1()
	cov0 = 1;
	scale = 1;
	p = 2;
	MPerDim = 5;
	Sigma1 = PowerExponentialCovarianceMatrix1D(cov0,scale,p,MPerDim);
	Sigma2 = PowerExponentialCovarianceMatrixAnyD(cov0,scale,p,MPerDim);
	if (Sigma1~=Sigma2)
		warning('Test1 failed');
	else
		disp('Test1 ok');
	end
end

function test2()
	cov0 = 1;
	scale = [1 2];
	p = 2;
	MPerDim = 5;
	Sigma1 = PowerExponentialCovarianceMatrix2D(cov0,scale,p,MPerDim);
	Sigma2 = PowerExponentialCovarianceMatrixAnyD(cov0,scale,p,MPerDim);
	if (Sigma1~=Sigma2)
		warning('Test2 failed');
	else
		disp('Test2 ok');
	end
end
