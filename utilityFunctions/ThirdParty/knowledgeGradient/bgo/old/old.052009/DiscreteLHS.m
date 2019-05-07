% design = DiscreteLHS(MPerDim,dim,n0)
% 
% Outputs n0 unique dim-dimensional points whose components are integers
% between 0 and MPerdim-1 inclusive.
function design = DiscreteLHS(MPerDim,dim,n0)

% Each output row needs to have each component an integer >= 0 and < MPerDim.
% The output of lhsdesign is points of the form (k+.5)/n0, where k is an
% integer.  By subtracting .5/n0, we get points of the form k/n0.
tries = 1;
design = helper(MPerDim,dim,n0);
while (size(design,1)<n0)
	if (tries > 20)
		warning(sprintf('Having trouble allocating a Latin Hypercube with no duplicates, tries=%d ', tries));
	end
	design = helper(MPerDim,dim,n0);
	tries = tries+1;
end
% Debugging
disp(sprintf('tries=%d', tries));
end

function design = helper(MPerDim,dim,n0)
	design = floor(MPerDim*(lhsdesign(n0,dim,'smooth','off')-.5/n0));
	design = unique(design,'rows');
end
