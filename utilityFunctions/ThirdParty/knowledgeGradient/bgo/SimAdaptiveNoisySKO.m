%[impl,x,y,extras] = SimAdaptiveNoisySKO(N,sample,dim,n0a,n0b,MPerDim,scale0,g0)
function [impl,x,y,extras] = SimAdaptiveNoisySKO(N,sample,dim,n0a,n0b,MPerDim,scale0,g0)
	scaleL = scale0./10;
	scaleU = scale0.*10;
	gL = 1e-10;
	gU = 1;
	MLE = @(xd,y)EstimateHyperparameters(xd,y,scale0,scaleL,scaleU,g0,gL,gU);
	[impl x y extras] = SimAdaptiveSKO(N,sample,dim,n0a,n0b,MPerDim,MLE);
end
