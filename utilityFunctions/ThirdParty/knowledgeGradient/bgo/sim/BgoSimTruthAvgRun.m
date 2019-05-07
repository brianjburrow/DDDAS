% Run with
% state=BgoSimTruthAvgRun(state,trueValList(state.nRuns+1,:))
function s=BgoSimTruthAvgRun(s,thisTrueVal)
    started = datestr(now());

    thisSample = @(x) double(thisTrueVal(x)) + normrnd(0,s.noise_sdev);
    best_truth = max(thisTrueVal);

    n = s.nRuns + 1;
    disp(sprintf('Run n=%d',n));
    tic;
    [impl,x,y,extras] = s.simulator(s.N,thisSample);
    %impl = s.simulator(s.N,thisSample);
    cpu = toc;
    assert(length(impl)==s.N)
    s.OC(n,:) = best_truth - thisTrueVal(impl);
    s.nRuns = s.nRuns+1;
    % save the extra stuff?
    s.CPU(n) = cpu;
    s.x(n,:) = x;
    s.y(n,:) = y;
    s.impl(n,:) = impl;
    s.extras{n} = extras;
   
  
    finished = datestr(now());
    s.modified = finished;
    disp(sprintf('Started at %s, Finished at %s', started, finished));
end
