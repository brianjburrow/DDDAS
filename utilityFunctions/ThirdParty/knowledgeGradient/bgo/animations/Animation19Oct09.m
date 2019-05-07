% This animation was created for a presentation to Cornell ORIE in September
% 2009.  The figures generated are 

% add subdirectories containing matlabKG libraries.
% addpath(genpath('~/work/cvs/matlabKG'));

% Place to save the figures generated.
figdir_name = '~/work/cvs/matlabKG/bgo/animations/Animation19Oct09/';

% check path
% more on; path; more off;

% Parameters you might want to change
scale = 30;
cov0 = 0.5;
noisevar = 0.2;
n_max = 50;
M=300;

s = RandStream.create('mt19937ar','seed',104);
RandStream.setDefaultStream(s);

belief_axis = [1 M -2 2.5];
kg_axis = [1 M -15 -1];

mu0 = 0;
p = 2;
noisestd = sqrt(noisevar);
Sigma0 = PowExpCov(cov0,scale,p,M,1);
y_true=mvnrnd(mu0*ones(1,M),Sigma0)';

x = [];
y = [];
n = 0;
while n < n_max 

  xdecision = PlotPosterior(y_true, cov0, scale, mu0, noisevar, x, y, belief_axis);

  figure(1);
  fig_name = [figdir_name,num2str(n),'.pdf'];
  print('-dpdf', fig_name )
  pause(1)
  % close(1);

  % Augment data

  x = [x; xdecision];
  y = [y; y_true(xdecision) + noisestd * randn()];
  n = length(x)
end
