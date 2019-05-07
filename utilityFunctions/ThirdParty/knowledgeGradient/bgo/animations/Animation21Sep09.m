% This animation was created for a presentation to Cornell ORIE in September
% 2009.  The figures generated are 

% add subdirectories containing matlabKG libraries.
addpath(genpath('~/work/cvs/matlabKG'));

% Place to save the figures generated.
figdir_name = '~/work/cvs/matlabKG/bgo/animations/Animation21Sep09/';

% check path
% more on; path; more off;

% Parameters you might want to change
scale = 30;
cov0 = 0.5;
noisevar = 0.2;
n_max = 20;
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

%  maximin LHD
 x = [50 150 250]';
 y = y_true(x) + noisestd * randn(size(x));

for n=0:length(x)
  n 
  subplot(2,1,1)
  [true_istar_val, true_istar] = max(y_true);
  plot(1:M, y_true, '-k', true_istar, true_istar_val, 'sk', x(1:n), y(1:n), 'bo');
  axis(belief_axis);
  xlabel('x');
  ylabel('value');

  subplot(2,1,2)
  axis(kg_axis)
  cla
  fig_name = [figdir_name,'pre.kg.',num2str(n),'.pdf'];
  print('-dpdf', fig_name )
  pause
end


n = length(x)
while n < n_max 

  xdecision = PlotKG(y_true, cov0, scale, mu0, noisevar, x, y, 0, belief_axis, kg_axis);

  figure(1);
  fig_name = [figdir_name,'kg.',num2str(n),'.pdf'];
  print('-dpdf', fig_name )
  pause
  % close(1);

  % Augment data

  x = [x; xdecision];
  y = [y; y_true(xdecision) + noisestd * randn()];
  n = length(x)
end
