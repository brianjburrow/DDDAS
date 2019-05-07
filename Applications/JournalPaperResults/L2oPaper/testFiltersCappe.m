clear all
close all
clc

%% testFiltersSimple
% This is the function used in "Nonlinear Kalman Filtering via Measure
% Change" to demonstrate the performance of each filter on a nonlinear
% state space model.
% numRepts controls the monte carlo analysis (i.e., how many times we
% repeat the simulation)
% idx = 2:someNumber controls the number of iterations we run each filter.


countFailUKF = 0;
countFailEKF = 0;
countFailPF = 0;
countFailPF2 = 0;
countFailL2O = 0;
numRepts = 1;
for dmx = 1:numRepts
% Initialize Problem
initialState = 0;
processNoise = 10;
measurementNoise = 1;
%% Initialize UKF

ukf = unscentedKalmanFilter( ...
@cappe_ukfStateFunc,...                                                     % Assign the state evolution function
@cappe_ukfDataGen, ...                                                      % Assign the forward model
initialState, ...                                                           % Set the initial conditions
'HasAdditiveMeasurementNoise',1);                                           % Let MATLAB know the type of error model we are using 
ukf.ProcessNoise = processNoise;                                            % Set the uncertainty of our evolution model
ukf.MeasurementNoise = measurementNoise;                                    % Set the uncertainty of our forward model
ukf.Alpha = 1;                                  
ukf.Kappa = 3;                                                              % These hyperparameters needed to be tuned to get the UKF to work well on this problem.
ukf.Beta = 8;                                                               % They control the "spread" of the sigma points

%% Initialize EKF
ekf = extendedKalmanFilter( ...                                             
@cappe_ukfStateFunc,...                                                     % Assign the state evolution function
@cappe_ukfDataGen, ...                                                      % Assign the forward model
initialState, ...                                                           % Set the initial conditions
'HasAdditiveMeasurementNoise',1);                                           % Let MATLAB know the type of error model we are using
ekf.ProcessNoise = processNoise;                                            % Set the evolution model uncertainty (a.k.a. process model uncertainty)
ekf.MeasurementNoise = measurementNoise;                                    % Set the forward model uncertainty


%% Initiazile Bootstrap Particle Filter
numParticles = 200;                                                         % Set the number of particles used at each iteration.  Increase for better accuracy but higher computational cost
% xx, yy, zz are x_t, x_t-1, time respectively
% xx, dat are current state and data respectively
pdf_g = @(dat,xx) normpdf(dat, 0.05 * xx^2, 1);                             % Set up the likelihood model

PF.resampling = 1;                                                          % Set resampling = true, this will change dynamically later
PF.likelihood = pdf_g; 
PF.proposalSampler.random = @(xx, dat, Tim) normrnd(...                     % Create a random number generator to sample from the state evolution function
    cappe_benchmark.StateTransitionC(xx, Tim),...
    processNoise^0.5);
PF.particles = zeros(length(initialState), numParticles);                   
PF.weights = zeros(1, numParticles);
for idx = 1:numParticles
    data = cappe_benchmark.GetData(initialState);                           %
    PF.particles(:,idx) = PF.proposalSampler.random(initialState, ...       % generate an initial particle
        data, 0);
    PF.weights(idx) = PF.likelihood(data, initialState);                    % Compute initial weights
end
PF.weights = PF.weights/sum(PF.weights);                                    % Normalize the weights



%% Initialize L2O_Filter

L2O.samples = linspace(-60, 60, 2000);                                      % Construct the offline library inputs
L2O.stateEvolution = @(xx, yy) cappe_benchmark.StateTransition(xx, yy);     % Set a random number generator for the state evolution
L2O.oper = 0;                                                               % Operating condition (time in this case)
L2O.output = cappe_benchmark.GetDataC(L2O.samples);                         % Compute offline forward model evaluations
L2O.R = measurementNoise;                                                   % Set the measurement noise
L2O.estimate = initialState;                                                % Set the initial state
L2O.err = processNoise;                                                     % Set the initial uncertainty
%% Run Problem

disp(dmx)
clear err err2 err3 err4 err5 
clear est est2 est3 est4 est5
clear state
    
state(1) = initialState;
est(1) = initialState;
err(1) = 3*10^0.5;
est2(1) = initialState;
err2(1) = 3*10^0.5;
est3(1) = initialState;
err3(1) = 3*10^0.5;
est4(1) = initialState;
err4(1) = 3*10^0.5;
est5(1) = initialState;
err5(1) = 3*10^0.5;


for idx = 2:70
    % Evolve State and collect data
    state(idx) = cappe_benchmark.StateTransition(state(idx - 1), idx);      % evolve the true state, this is what we are trying to predict accurately.
    data = cappe_benchmark.GetData(state(idx));                             % Collect a noisy measurement from the system
    % Evaluate UKF
    tic;
    [est(idx), err(idx)] = correct(ukf, data, idx);
    predict(ukf, idx);
    ukfTime(idx) = toc;
    
    % Evaluate the Bootstrap Particle Filter
    PF.data = data;
    PF.oper = idx;
    tic
    PF = BootstrapFilter(PF);
    pfTime(idx) = toc;
    
    est2(idx) = PF.estimate;
    err2(idx) = PF.err;
   
    
    %Evaluate L2O_filter
    L2O.oper = idx;
    L2O.data = data;
    tic
    L2O = L2oFilter_up(L2O, 0, 50, 500, 1);
    l2oTime(idx) = toc;
    est3(idx) = L2O.estimate;
    err3(idx) = L2O.err;
    disp('time')
    disp(l2oTime(idx))
    
    
    % Evaluate EKF
    tic
    [est4(idx), err4(idx)] = correct(ekf, data, idx);
    predict(ekf, idx);
    ekfTime(idx) = toc;
    
    if abs(est(idx) - state(idx)) > 3*err(idx)^0.5
        countFailUKF = countFailUKF + 1;
    end
    if abs(est2(idx) - state(idx)) > 3*err2(idx)^0.5
        countFailPF = countFailPF + 1;
    end
    if abs(est3(idx) - state(idx)) > 3*err3(idx)^0.5
        countFailL2O = countFailL2O + 1;
    end
    if abs(est4(idx) - state(idx)) > 3*err4(idx)^0.5
        countFailEKF = countFailEKF + 1;
    end

end

for tmz = 1:70
    E1(dmx, tmz) = sqrt(mean((est(:,tmz)   -  state(:,tmz)).^2));
    E2(dmx, tmz) = sqrt(mean((est2(:,tmz)  -  state(:,tmz)).^2));
    E3(dmx, tmz) = sqrt(mean((est3(:,tmz)  -  state(:,tmz)).^2));
    E4(dmx, tmz) = sqrt(mean((est4(:,tmz)  -  state(:,tmz)).^2));
    S1(dmx, tmz) = mean(mean(err));
    S2(dmx, tmz) = mean(mean(err2));
    S3(dmx, tmz) = mean(mean(err3));
    S4(dmx, tmz) = mean(mean(err4));
end

end

EE1 = mean(mean(E1'));
EE2 = mean(mean(E2'));
EE3 = mean(mean(E3'));
EE4 = mean(mean(E4'));

EEE1 = mean(mean(E1));
EEE2 = mean(mean(E2));
EEE3 = mean(mean(E3));
EEE4 = mean(mean(E4));

h = figure;
set(gcf,'color','w');
subplot(2, 2, 1);
scatter(linspace(1,length(est),length(est)), est-state, 'black', 'filled','LineWidth', 0.05)
hold on
errorbar(zeros(1,length(est)), 3*err.^0.5, 'blue','LineStyle',  '-')
xlabel('Time')
ylabel('Error')
ylim([-60,60])
xlim([0, 80])
title('UKF')

subplot(2, 2, 2);
scatter(linspace(1,length(est2),length(est2)), est2-state, 'black', 'filled','LineWidth', 0.05)
hold on
errorbar(zeros(1,length(est2)), 3*err2.^0.5, 'blue','LineStyle',  '-')
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 80])
title('Particle Filter')

subplot(2, 2, 3);
scatter(linspace(1,length(est3),length(est3)), est3-state, 'black','filled', 'LineWidth', 0.05)
hold on
errorbar(zeros(1,length(est3)), 3*err3.^0.5, 'blue','LineStyle',  '-')
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 80])
title('L2O Filter')

subplot(2, 2, 4);
scatter(linspace(1,length(est4),length(est4)), est4-state, 'black', 'filled', 'LineWidth', 0.05)
hold on
errorbar(zeros(1,length(est4)), 3*err4.^0.5, 'blue','LineStyle',  '-')
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 80])
title('EKF')


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

h = figure;
set(gcf,'color','w');
subplot(2, 2, 1);
plot(est, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
errorbar(est, 3*err.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est),length(est)), state, 'black', 'filled', 'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('State')
ylim([-40, 50])
xlim([0, 70])
title('UKF')

subplot(2, 2, 2);
errorbar(est2, 3*err2.^0.5, 'blue','LineStyle',  '-')
hold on
plot(est2, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
scatter(linspace(1,length(est2),length(est2)), state, 'black', 'filled', 'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('State')
ylim([-40, 50])
xlim([0, 70])
title('Particle Filter')

subplot(2, 2, 3);
errorbar(est3, 3*err3.^0.5, 'blue','LineStyle',  '-')
hold on
plot(est3, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
scatter(linspace(1,length(est3),length(est3)), state, 'black','filled', 'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('State')
ylim([-40, 50])
xlim([0, 70])
title('L2O Filter')
hold on

subplot(2, 2, 4);
plot(est4, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
scatter(linspace(1,length(est4),length(est4)), state, 'black','filled',  'LineWidth', 0.05)
hold on
errorbar(est4, 3*err4.^0.5, 'blue','LineStyle',  '-')
xlabel('Time')
ylabel('State')
ylim([-40, 50])
xlim([0, 70])
title('EKF')
hold on

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
hold on