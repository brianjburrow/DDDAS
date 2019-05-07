function [w, feval] = frankwolfe(X,b, w, Hw, feval, itr) 
  n = size(X,1);
  d = size(X,2);

  H = zeros(n,1);
  w0 = zeros(n,1);
  g1 = Hw;
  g  = g1 - b;
  for i=1:n
    feval(1) = feval(1) + w(i)*(0.5*g1(i)-b(i));
  end

  ii = 1;
  while (ii<itr+1)
   
    ii = ii+1;
    [~,I] = min(g);
    w0(I) = 1;
    
    c1 = 0;
    c4 = 0;

  % Build column of H(:,I)
    for j=1:n  
      H(j) = 1.0;
      for k=1:d  
          if X(I,k) < X(j,k)
                H(j) = H(j)*(1-X(j,k));
          else
                H(j) = H(j)*(1-X(I,k));
          end
      end
      c1 = c1 + 0.5*w(j)*g1(j);
      c4 = c4 + w(j)*b(j);
    end
    
    c2 = g1(I);
    c3 = H(I)/2;
    c5 = b(I); 
    p1 = (c1-c2+c3);
    p2 = (c4-2*c1-c5+c2);
    p3 = (c1-c4);
    
    if 1 < -p2/2/p1  
      G = 1;
    else  
      %G = -p2/2/p1;
      G = 2/(2 + ii);
    end
    
    for j=1:n 
      w(j)  = w(j) + G*(w0(j)-w(j));
      g1(j) = (1-G)*g1(j) + G*H(j);
      g(j)  = g1(j) - b(j);
    end
    
    for j=1:n  
      feval(ii) = feval(ii) + w(j)*(0.5*g1(j)-b(j));
    end

    w0(I) = 0;
  end
end
