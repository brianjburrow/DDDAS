% state = BgoSimRun(old_state, nRuns)
% 
% Runs the simulation specified by the parameters in the passed state for the
% specified number of runs, and the state return has the new results appended
% to the old results in the passed state.  The state should be created by
% BgoSimCreate.m.
function s = BgoSimRun(old_s, nRuns)
    started = datestr(now());
    nRunsAlready = old_s.nRuns;

    s = old_s;
    s.CPU = [old_s.CPU; zeros(nRuns,1,'single')];
    s.OC = [old_s.OC; zeros(nRuns,s.N)];
    s.x = [old_s.x; zeros(nRuns,s.N, 'uint16')];
    s.y = [old_s.y; zeros(nRuns,s.N, 'single')];
    s.impl = [old_s.impl; zeros(nRuns,s.N, 'uint16')];
    s.extras = old_s.extras; % Can I preallocate a cell array?

    for n=1+nRunsAlready:nRunsAlready+nRuns
	disp(sprintf('Run n=%d',n));
	tic;
	[impl,x,y,extras] = s.simulator(s.N,s.sample);
	s.CPU(n) = toc;
	for k=1:length(impl)
		s.OC(n,k) = s.best_truth - s.truth(impl(k));
	end
	s.x(n,:) = x;
	s.y(n,:) = y;
	s.impl(n,:) = impl;
	s.extras{n} = extras;
    end
    s.nRuns = n;
  
    finished = datestr(now());
    s.modified = finished;
    disp(sprintf('Started at %s, Finished at %s', started, finished));
end
