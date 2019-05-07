function TestBayesUpdateBlock()
test1();
test2();
end

function test1()
MPerDim=6;
dim = 2;
scale = [1 1];
cov0 = 1;
M = MPerDim^dim;
p = 2;
Sigma0=PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
mu0 = zeros(M,1);
truth = rand(M,1);
x = [1:4:M]';
y = truth(x);
noises = zeros(size(x));
[mu,Sigma]=BayesUpdateBlock(mu0,Sigma0,x,y,noises);
[mu1,Sigma1]=helper1(mu0,Sigma0,x,y,noises);
[mu2,Sigma2]=helper2(mu0,Sigma0,x,y,noises);
[mu3,Sigma3]=BayesUpdateBlock(mu0,Sigma0,repmat(x,2,1),repmat(y,2,1),repmat(noises,2,1));
[mu3,Sigma3]=BayesUpdateBlockNaive(mu0,Sigma0,repmat(x,2,1),repmat(y,2,1),repmat(noises,2,1));
mudiff1 = full(max(abs(mu1-mu)));
Sigmadiff1 = max(max(abs(Sigma1-Sigma)));
mudiff2 = full(max(abs(mu2-mu)));
Sigmadiff2 = max(max(abs(Sigma2-Sigma)));
mudiff3 = full(max(abs(mu3-mu)));
Sigmadiff3 = max(max(abs(Sigma3-Sigma)));
if (max([mudiff1,mudiff2,Sigmadiff1,Sigmadiff2,mudiff3,Sigmadiff3])>1e-10)
    warning('The difference between different methods of updating mu and Sigma is large.  There may be a bug.');
	disp([mudiff1, mudiff2, mudiff3]);
	disp(full([Sigmadiff1, Sigmadiff2, Sigmadiff3]));
	(mu-mu1)'
	(mu-truth)'
	(mu1-truth)'
else
    disp('Test1 ok');
end
end

function test2()
MPerDim=6;
dim = 2;
scale = [1 1];
cov0 = 1;
M = MPerDim^dim;
p = 2;
Sigma0=PowerExponentialCovarianceMatrix(cov0,scale,p,MPerDim);
mu0 = zeros(M,1);
truth = rand(M,1);
x = 1 + floor(M*rand(1,20))';
y = truth(x) + normrnd(0,1);
noises = ones(size(x));
[mu,Sigma]=BayesUpdateBlock(mu0,Sigma0,x,y,noises);
[mu1,Sigma1]=helper1(mu0,Sigma0,x,y,noises);
[mu2,Sigma2]=helper2(mu0,Sigma0,x,y,noises);
mudiff2 = max(abs(mu2-mu));
Sigmadiff2 = max(max(abs(Sigma2-Sigma)));
mudiff1 = max(abs(mu1-mu));
Sigmadiff1 = max(max(abs(Sigma1-Sigma)));
if (max([mudiff1,mudiff2,Sigmadiff1,Sigmadiff2])>1e-10)
    warning('The difference between different methods of updating mu and Sigma is large.  There may be a bug.');
else
    disp('Test2 ok');
end
end


function [newmu,newSigma]=helper1(mu,Sigma,x,y,noises);
	newmu = mu;
	newSigma = Sigma;
	for i=1:length(x)
		[newmu,newSigma]=BayesUpdate(newmu,newSigma,x(i),y(i),noises(i));
	end
end

function [newmu,newSigma]=helper2(mu,Sigma,x,y,noises);
	newmu = mu;
	newSigma = Sigma;
	% Note that this will produce error messages.
	for i=1:length(x)
		[newmu,newSigma]=BayesUpdateBlock(newmu,newSigma,x(i),y(i),noises(i));
	end
end
