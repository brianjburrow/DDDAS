function [sample, accept, C] = gpMetropolisHastings(currentState, proposalState,...
    propPdf, targetPdf, iteration, C0, samples)

% All probability distributions passed must return log probabilities
% propPdf is the proposal probability density function.  Pass as a function
%            handle.  Should take 2 inputs, the current state and proposal
%            . state
% Likewise for the targetPdf, but only one input

%% Computing the acceptance ratio
t = iteration;
t0 = 5000;
sd = (2.4^2)/length(currentState);
epsilon = 10^-14;

[muProp, stdProp] = targetPdf(proposalState);
a = normrnd(muProp, stdProp, 1);
b = propPdf(currentState, proposalState, C0);
[muCur, stdCur] = targetPdf(currentState);
c = normrnd(muCur, stdCur, 1);
d = propPdf(proposalState, currentState, C0);

log_a = a + b - c - d;


alpha = min(log(1), log_a);


alpha = exp(alpha);
u = rand();


if u <= alpha
    sample = proposalState;
    accept = true;
else
    sample = currentState;
    accept = false;
end



if iteration <= t0
    C = C0;
else
    C = sd*cov(samples) + sd*epsilon*eye(length(C0));
end


end