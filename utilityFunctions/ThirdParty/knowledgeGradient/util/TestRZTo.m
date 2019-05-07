%{
Unit test for the *To*.m functions.  If there is an error, an assertion will
fail.  Otherwise it will print out a 'Test OK' message.
%}
delta = .01;
M=100;
d=5;
origin = ones(1,5);
x=[1.1,1.0,1.98,1.05,1.99];
zd=RdToZd(x,delta,origin);
z=ZdToZ(zd,M,d);
zd2=ZToZd(z,M,d);
assert(all(zd2==zd));
x2=ZdToRd(zd,delta,origin);
assert(all(x2==x));
display('Test OK');
