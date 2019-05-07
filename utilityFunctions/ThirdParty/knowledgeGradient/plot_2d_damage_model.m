close all
clc
%% Description
% This code is used to perform offline optimization of a vehicle capability
% model for a 2-D toy physics problem.  The capability model is described 
% as a Gaussian process model from strains to flutter speed.  
 
% The damage problem takes 2-D physics inputs: 
% --- damage takes values between 0 and 1, and represents normalized crack
% -density.
% --- loadFactor takes values between 1.0 and 6.0, and represents a uniformly
% -distributed load on a square cantilever beam, which is solved via Euler
% -Bernoulli.  See missionControl.py for details on the damage and bea
% -models.

% The optimization problem: Goal is to select the training cases used to
% build the GP capability model in an optimal way, assuming some limit on
% the number of training data that can be used online, subject to some
% prior information about the likely damage/loadFactor pairings.

% The optimization problem's objective is to maximize the information gain:
% described by the K-L divergence between the posterior and the prior at
% each stage.  This is computed via Monte Carlo integration over the
% current
%% User Input Block
% Set up damage problem
start_number = 30;
damages = linspace(0, 1, start_number);  %  training inputs
loadFactor = linspace(1, 6, start_number); % Can vary this as well to make this a 2-d problem
gauge_loc = 24; % meters
noise_level = 0.001; % strain gauge noise
num_MC_repeats = 50; % collect M noisy strains, calc entropy, repeat num_MC_repeats to get avg entropy
num_noisy_strains_per_damage = 100;
prior_std = 0.1;
% Set up KL divergence problem
xx = linspace(0, 1, 10); % set up training points
xxStrains = zeros(1,length(damages));
xxFlutterSpeeds = zeros(1,length(damages));

% Generate original training samples for the capability model
strains = zeros(length(loadFactor),length(damages));
flutterSpeeds = zeros(length(loadFactor),length(damages));
for dmx = 1:start_number
    for idx = 1:start_number
        strains(idx, dmx) = py.mymod2.pass_strains_to_matlab(...
            damages(idx), loadFactor(dmx), gauge_loc...
            ); %  noiseless training outputs
        flutterSpeeds(idx, dmx) = py.mymod2.pass_flutter_to_matlab(...
            damages(idx), loadFactor(dmx), gauge_loc...
            );
    end
end

training_inputs = zeros(2, start_number^2);
training_outputs = zeros(1,start_number^2);

counter = 1;
for idx = 1:start_number
    for dmx = 1:start_number
        training_inputs(1, counter) = loadFactor(idx);
        training_inputs(2, counter) = strains(dmx, idx);
        training_outputs(counter) = flutterSpeeds(dmx, idx);
        counter = counter + 1;
    end
end

subplot(3, 1, 1)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Times','fontsize',12)
[X, Y] = meshgrid(loadFactor, damages);
Z = strains;
disp(strains)

contourf(X, Y, Z, 10);
xlabel('Load Factor', 'fontsize', 20)
ylabel({'Normed' ; 'Crack Density'}, 'fontsize', 20)
zlabel('Strain', 'fontsize', 20)
colorbar;

legend('Strain','Location','Best')

subplot(3, 1, 2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Times','fontsize',12)

[X, Y] = meshgrid(loadFactor, damages);
Z = flutterSpeeds;

contourf(X, Y, Z, 10);

xlabel('Load Factor', 'fontsize', 20)
ylabel({'Normed' ; 'Crack Density'}, 'fontsize', 20)
zlabel('Flutter Speed', 'fontsize', 20)
set(gcf,'color','w');
legend('Flutter Speed','Location','Best')
colorbar;


subplot(3, 1, 3)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Times','fontsize',12)
% make up some data
x = training_inputs(2,:);
y = training_inputs(1,:);
z = training_outputs;
n = 1000;
[X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));

%create contour plot
contourf(X,Y,griddata(x,y,z,X,Y), 10)


ylabel('Load Factor', 'fontsize', 20)
xlabel('Strains', 'fontsize', 20)
zlabel('Flutter Speed', 'fontsize', 20)
set(gcf,'color','w');
legend('Flutter Speed','Location','Best')
colorbar;








