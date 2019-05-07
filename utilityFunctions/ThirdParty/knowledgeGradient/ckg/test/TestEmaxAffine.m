% Test code for EmaxAffine
function ok = TestEmaxAffine()
ok = 1;
ok = test1() && ok;
ok = test2() && ok;
ok = test3() && ok;
ok = test4() && ok;
ok = test5() && ok;
ok = test6() && ok;
for i=[1:100], test7() && ok; end
ok = 1; for i=[1:100], ok = test8() && ok; end
if (ok), disp('test8: OK'); else, disp('test8: ERR'); end
test9();
end

% Calculates for zero input and a single nonzero slope of 1.
% The answer should be the same as Exp Z^+ = 1/sqrt(2Pi) = 
function ok = test1()
M=5;
a=zeros(M,1);
b=zeros(M,1);
b(1)=1;
result = EmaxAffine(a,b);
should = 1/sqrt(2*pi);
display(sprintf('TestEmaxAffine: test1: result=%f should=%f', result, should));
ok = abs(result-should)<1e-8;
end

function ok = test2()
b = rand(5,1);
a = zeros(5,1);
result = EmaxAffine(a,b);
should = (max(b)-min(b))/sqrt(2*pi);
display(sprintf('test2: result=%f should=%f', result, should));
ok = abs(result-should)<1e-8;
end

function ok = test3()
a = [ -0.6712,   -0.2684,    0.1389,    0.2616,    0.0531,   -0.0500,    0.1266,    0.3320,    0.5411,    0.7563]';
b = [  0.0053,    0.0505,    0.1070,    0.1872,    0.1168,    0.0723,    0.0438,    0.0344,    0.0253,    0.0165]';
ok = MCresult('test3',a,b);
end

function ok = test4()
% This test is based on the fact that, in general, 
% EmaxAffine([0;0],[0;const]) = \Exp[(const*Z)^+] = |const|/sqrt(2*pi).
const = 5;
result = EmaxAffine([0;0],[0;const]);
should = abs(const)/sqrt(2*pi);
display(sprintf('test4: result=%f should=%f', result, should));
ok = abs(result-should)<1e-8;
end

function ok = test5()
a = [  -0.6712,   -0.2684,    0.1389,    0.2616,    0.0531,   -0.0500,    0.1266,    0.3320,    0.5411,    0.7563]';
b = [   0.0234,    0.0146,    0.0092,    0.0057,    0.0036,    0.0022,    0.0013,    0.0011,    0.0008,    0.0005]';
ok = MCresult('test5',a,b);
end

function ok = test6()
a = [  -0.6712,   -0.2684,    0.1389,    0.2616,    0.0531,   -0.0500,    0.1266,    0.3320,    0.5411,    0.7563]';
b = [   0.0009,    0.0090,    0.0190,    0.0332,    0.0548,    0.0886,    0.1420,    0.2570,    0.1891,    0.1233]';
ok = MCresult('test6',a,b);
end

% This test is different in character than the other ones, because it
% generates an a and b at random, and then compares to a monte-carlo
% evaluation of the expectation.
function ok = test7()
nsamples = 10^4;
M=ceil(10*rand);
a = rand(10,1);
b = rand(10,1);
ok = MCresult('test7',a,b);
end


function ok = MCresult(testname, a, b)
    result = EmaxAffine(a,b);
    [should,err] = EmaxAffineMonteCarlo(a,b,10^4);
    if (result > should - 3*err && result < should + 3*err) 
        display(sprintf('%s: result=%f should=%f+/-%f OK', testname, result, should, err));
	ok = 1;
    else
        display(sprintf('%s: result=%f should=%f+/-%f ERR', testname, result, should, err));
        % Rerun it with more samples, in case it was just a monte carlo cosmic ray.
        [should,err] = EmaxAffine(a,b,10^5);
        if (result > should - 3*err && result < should + 3*err) 
            display(sprintf('rerun: should=%f+/-%f OK', should, err));
	    ok = 1;
        else
            % Display the inputs so that we can investigate the case.
            display(sprintf('rerun: should=%f+/-%f ERR', should, err));
            a=a'; b=b';
            display(a);
            display(b);
	    ok = 0;
        end
 
    end
end

% Generates an a and b at random, and then compares the current version of
% EmaxAffineBreakpoints with a previous version.
function ok=test8()
M=ceil(10*rand);
a = rand(10,1);
b = rand(10,1);
[c1, keep1] = AffineBreakpoints(a,b);
[c2, keep2] = AffineBreakpointsV2(a,b);
if (length(c1)~=length(c2) || any(c1 ~= c2) || length(keep1) ~= length(keep2) || any(keep1~=keep2))
	disp(sprintf('test8: FAIL a=%s b=%s',mat2str(a),mat2str(b)));
	ok = 0;
else
	ok = 1;
	%disp('test8: OK');
end
end


function ok=test9()
M=10^5;
a = rand(M,1);
b = rand(M,1);
tic, result = EmaxAffine(a,b);
t = toc;
disp(sprintf('test9: (%d seconds for M=%d) OK', t, M));
end
