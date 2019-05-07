% function [mu,sposteriorEGO,snoise] = TestEGOSKO(n)
% Tests the EGO and SKO functions.  Returns some debugging information on
% failure.  The argument n is the number of trials to run, and defaults to 100
% if you don't specify it.
function [mu,sposteriorEGO,snoise] = TestEGOSKO(n)
    if (nargin==0)
	n=100;
    end
    for i=1:n
	[mu,sposteriorEGO,sposteriorSKO,snoise,xEGO,xSKO]=TestEGOSKOHelper();
	if (xEGO~=xSKO)
		disp('Test failed');
		% Display for debugging.
		mu
		sposteriorEGO
		snoise
		xEGO
		xSKO
		return
	end
    end
    disp('Test OK');
end



function [mu,sposteriorEGO,sposteriorSKO,snoise,xEGO,xSKO] = TestEGOSKOHelper()
% Parameters to tweak
M=100;
percent_measured = .2;
snoise = .015; % Should be a strictly positive number close to 0

sposterior = rand(1,M);
% If an unmeasured sposterior is close to snoise, SKO and EGO need not
% necessarily agree. We avoid that contingency here.
sposterior(find(sposterior<10*snoise))=.5;
mu = rand(1,M);
measured = find(rand(1,M)<percent_measured);
% EGO is not designed to handle the case when we have measured nothing or
% everything.  SKO is not designed to handle the case when we have measured
% nothing.  The code works for those cases, but spits out a warning.  We skip
% over that warning by doing the following.
while (length(measured)==0 || length(measured)==M)
	measured = find(rand(1,M)<percent_measured);
end
sposteriorEGO = sposterior;
sposteriorSKO = sposterior;
sposteriorEGO(measured) = 0;
sposteriorSKO(measured) = snoise;
xEGO = EGO(mu', sposteriorEGO');
xSKO = SKO(mu', sposteriorSKO', snoise, measured);
disp(sprintf('xEGO=%d xSKO=%d', xEGO, xSKO));
end
