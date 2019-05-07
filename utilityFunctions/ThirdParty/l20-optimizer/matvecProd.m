function result = matvecProd(m1, v1)
n = size(m1,1);
d = size(m1,2);

  result = zeros(n,1);
   
  for ii=1:n
    for jj=1:n
      A = 1.0;
      for kk=1:d
          if (m1(jj,kk) < m1(ii,kk)) 
            A = A* (1-m1(ii,kk)) ;
          else
            A = A* (1-m1(jj,kk));
          end
      end
      result(ii,1) = result(ii,1) + v1(jj)*A;
    end
  end

end