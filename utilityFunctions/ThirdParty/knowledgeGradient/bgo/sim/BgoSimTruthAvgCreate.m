% state = BgoSimTruthAvgCreate(N, M, truth, noise_sdev, simulator)
%
% This function creates a structure that contains the parameters of an
% experiment.  It can be passed to BgoSimTruthAvgRun(), which will actually run
% the experiment.  Compared to BgoSimCreate, this function creates a state for
% running simulations across multiple truths, while BgoSimCreate uses just one
% truth.   The truths to use are passed in one at a time to BgoSimTruthAvgRun.  We
% do not store them here because they are often really big, and we use them for
% multiple states, so we only want to store them once.
%
% simulator should be a function handle of the form
% implementation_decisions = simulator(N,sample);
 
function state = BgoSimTruthAvgCreate(N, M, noise_sdev, simulator, keepextras)
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
    state.noise_sdev = noise_sdev;
    state.simulator = simulator;
end
