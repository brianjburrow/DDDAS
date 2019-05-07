% state = BgoSimCreate(N, M, truth, noise_sdev, simulator)
%
% This function creates a structure that contains the parameters of an
% experiment.  It can be passed to BgoSimRun.m, which will actually run the
% experiment subject to various stopping conditions.
%
% simulator should be a function handle of the form
% implementation_decisions = simulator(N,sample);
 
function state = BgoSimCreate(N, M, truth, noise_sdev, simulator, keepextras)
    state.created = datestr(now());
    state.CPU = []; % CPU runtime
    state.OC = []; % Opportunity costs
    state.x = []; % Measurement decisions
    state.y = []; % Observations
    state.impl = []; % Implementation decisions
    if (nargin < 6)
	keepextras = 1;
    end
    state.keepextras = keepextras;
    if (keepextras)
	% extras may store adaptively estimate hyperparameters, or log-KG factors,
	% or anything else particular to the algorithm.  They are sometimes
	% large though, so storing them is optional.
	state.extras = {};
    end
    state.nRuns = 0;
    state.N = N;
    state.M = M;
    state.truth = truth;
    for i=1:M
	tmp(i) = truth(i);
    end
    state.best_truth = max(tmp);
    state.noise_sdev = noise_sdev;
    state.sample = @(x) truth(x) + normrnd(0,noise_sdev);
    state.simulator = simulator;
end
