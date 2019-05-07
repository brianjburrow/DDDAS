function [m1,m2] = unitHyper(m1,m2)  
  m1 = m1 - min([m1;m2]);
  m2 = m2 - min([m1;m2]);

  scale = max([m1;m2]);
  
  m1 = m1./scale;
  m2 = m2./scale;
end