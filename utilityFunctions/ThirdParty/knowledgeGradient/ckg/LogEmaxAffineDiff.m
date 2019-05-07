% This function calculates
% y = EmaxAffine(a1,b1)-EmaxAffine(a2,b2) in a numerically precise way.
% returned are logdiff = log(abs(y)) and sgndiff=sgn(y).

function [logdiff,sgndiff]=LogEmaxAffineDiff(a1,b1,a2,b2,c1,c2)
  a = ConvertToColumnVector(a);
  b = ConvertToColumnVector(b);
  if (length(a) ~= length(b))
	error('LogEmaxAffineDiff: a and b must be vectors of the same length');
  end

  if (nargin<6)
	[a1,b1] = AffineBreakpointsPrep(a1,b1);
	[c1, keep] = AffineBreakpoints(a1,b1);
	a1 = a1(keep);
	b1 = b1(keep);
	c1 = c1([1,keep+1]);
	[a2,b2] = AffineBreakpointsPrep(a2,b2);
	[c2, keep] = AffineBreakpoints(a2,b2);
	a2 = a2(keep);
	b2 = b2(keep);
	c2 = c2([1,keep+1]);
  end
		
  % Weave the c1 and c2 together to get a unified array.
  i1 = 1;
  i2 = 1;
  j = 1;

  cweave(1)=-Inf;
  while (i1<=length(a1) && i2<=length(a2))
  %a1weave(j)=a1(i1);
  b1weave(j)=b1(i1);
  %a2weave(j)=a2(i2);
  b2weave(j)=b2(i2);

  if (c1(i1+1)<c2(i2+1))
	cweave(j+1)=c1(i1+1);
	i1 = i1+1;
  elseif (c1(i1+1)>c2(i2+1))
	cweave(j+1)=c2(i2+1);
	i2 = i2+1;
  else
	cweave(j+1)=c1(i1+1);
	i1 = i1+1;
	i2 = i2+1;
  end
  j = j+1;
  end
  cweave(j)=Inf;
  assert(c1(i1)==Inf);
  assert(c2(i2)==Inf);
	
  %{
  disp(sprintf('c1=%s c2=%s', mat2str(c1),mat2str(c2)));
  disp(sprintf('b1=%s b2=%s', mat2str(b1),mat2str(b2)));
  disp(sprintf('cweave=%s', mat2str(cweave)));
  disp(sprintf('b1weave=%s b2weave=%s', mat2str(b1weave),mat2str(b2weave)));
  disp(sprintf('b1weave-b2weave=%s', mat2str(b1weave-b2weave)));
  %}
  %a = a1weave-a2weave
  b = b1weave-b2weave;
  c = cweave;

  bdiff = diff(b1weave-b2weave);
  if (all(bdiff==0))
	sgndiff=0;
	logdiff=-Inf;
  end
  logbdiff=log(abs(bdiff));
  pos=find(bdiff>0);
  neg=find(bdiff<0);

  if (length(pos)>0)
    	logsumpos = LogSumExp(logbdiff(pos)+LogEI(-abs(c(pos+1))));
  else
        logsumpos = -Inf;
  end
  if (length(neg)>0)
        logsumneg = LogSumExp(logbdiff(neg)+LogEI(-abs(c(neg+1))));
  else
        logsumneg = -Inf;
  end
    
  if (logsumpos>logsumneg)
	sgndiff=1;
	logdiff=LogMinusExp(logsumpos,logsumneg);
  elseif (logsumpos<logsumneg)
	sgndiff=-1;
	logdiff=LogMinusExp(logsumneg,logsumpos);
  else
	sgndiff=0;
	logdiff=-Inf;
  end

    %{
    logei1=LogSumExp(log(diff(b1weave))+LogEI(-abs(c(2:end-1)')));
    logei2=LogSumExp(log(diff(b2weave))+LogEI(-abs(c(2:end-1)')));
    d = log(abs(exp(logei1)-exp(logei2)));
    disp(sprintf('logei1=%g logei2=%g log(abs(exp(logei1)-exp(logei2)))=%g logdiff=%g sgndiff=%d', ...
	logei1, logei2, d, logdiff, sgndiff));
    %}
end
