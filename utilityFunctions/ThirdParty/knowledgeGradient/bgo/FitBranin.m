% Fits Branin 
function [cov0,theta] = FitBranin(MPerDim)
xmin = [-1.6,-.8];
xmax = [3.2, 1.6];
delta = (xmax-xmin)/MPerDim;
dim = 2;
toRd = @(z) ZdToRd(ZToZd(z,MPerDim,dim),delta,xmin);
for i=1:MPerDim^dim
	xd(i,:) = toRd(i);
	y(i)=Branin(xd(i,1),xd(i,2));
end
theta0 = [10 10]; lob = [1e-2 1e-2]; upb = [20 20];
[cov0,theta]=EstimateHyperparametersNoNoise(xd,y,@corrgauss,theta0,lob,upb);
end
