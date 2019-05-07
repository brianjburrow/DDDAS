% //**********************************************************************
% // Build H(R1:R2,:) and perform matrix-vector product
% //**********************************************************************
function result = buildXYb(m1, m2)
  %vector<double> partial_b(r2-r1,0.0);
   n = size(m1,1);
   m = size(m2,1);
   d = size(m1,2);
   result = zeros(n,1);
    for iI=1:n
        for jJ = 1:m
          A = 1.0;
          for kk=1:d
              if (m1(iI,kk) < m2(jJ,kk))
                    A = A* (1-m2(jJ,kk));
              else
                  A = A* (1-m1(iI,kk));
              end
          end
          result(iI) = result(iI)+A/m; 
        end
    end
end