function ok=LogEmaxAffineTest()
  ok = 1;
  for i=1:100
	ok = test1() && ok;
  end
  if (~ok)
    disp('LogEmaxAffineTest: FAIL');
  else 
    disp('LogEmaxAffineTest: OK');
  end
end

function ok = test1()
  M = 5;
  a = rand(M,1);
  b = rand(M,1);
  y1 = exp(LogEmaxAffine(a,b)) + max(a);
  y2 = EmaxAffine(a,b);
  ok = 1;
  if (abs(y1-y2)>10*eps(y2))
	ok = 0;
	disp(sprintf('error: a=%s b=%s y1=%f y2=%f', mat2str(a), mat2str(b), y1, y2));
  end
end
