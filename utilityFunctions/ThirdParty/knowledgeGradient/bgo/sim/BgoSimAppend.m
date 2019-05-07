function snew=BgoSimAppend(s,t)
	snew = s;
	snew.CPU = [s.CPU;t.CPU];
	snew.OC = [s.OC;t.OC];
	snew.x = [s.x;t.x];
	snew.y = [s.y;t.y];
	snew.impl = [s.impl;t.impl];
	k = length(s.extras);
	for n=1:length(t.extras)
		snew.extras{k+n} = t.extras{n};
	end
	snew.nRuns = s.nRuns + t.nRuns;
	snew.modified = datestr(now());
