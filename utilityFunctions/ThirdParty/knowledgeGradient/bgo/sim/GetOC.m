function [oc,err]=GetOC(OCarray,N,blocksize)
	data = MakeBlocks(OCarray(:,N),blocksize);
	nblocks = length(data)
	oc=mean(data);
	err=std(data)/sqrt(nblocks-1);
end
