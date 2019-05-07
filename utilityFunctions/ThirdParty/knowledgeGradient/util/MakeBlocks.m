% Takes a matrix input, and averages together blocks of rows.  So if "in" is an
% 1000x20 matrix, then MakeBlocks(in, 10) would return a 100x20 matrix "out",
% where out(i,j) would be the average of 10 entries from in(:,j). 
%
% This function is intended to be used for cases where each row of the input is
% the result of 1 run of a monte carlo simulation, and these MC results are not
% normally distributed.  By averaging the results together in blocks that are
% big enough, one hopes to obtain approximately normally distributed results.
% To test the normality of the block-averaged outputs in column N, one could
% run
% 	out = MakeBlocks(in, blocksize);
%	hist(out(:,N))
% and then look to see that the histogram is approximately normal.
function out = MakeBlocks(in, blocksize)
	if (blocksize <= 1)
		error('Must pass a blocksize > 1');
	end
	[nrows, ncols] = size(in);
	leftover = mod(nrows, blocksize);
	if (leftover ~= 0)
		disp(sprintf('Discarding %d rows', leftover));
	end
	nblocks = floor(nrows/blocksize);
	out = zeros(nblocks,ncols);
	for i=1:nblocks
		first = 1+(i-1)*blocksize;
		last = first + blocksize - 1;
		out(i,:) = mean(in(first:last,:));
	end
end
