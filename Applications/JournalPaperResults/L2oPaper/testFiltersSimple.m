clear all
clc
%% testFiltersSimple
% This is the function used in "Nonlinear Kalman Filtering via Measure
% Change" to demonstrate the performance of each filter on a linear state
% space model.
% numRepts controls the monte carlo analysis (i.e., how many times we
% repeat the simulation)
% idx = 2:someNumber controls the number of iterations we run each filter.
countFailUKF = 0;                                                           % Store the number of times the UKF's true error lies outside the error bounds.
countFailEKF = 0;
countFailPF  = 0;
countFailPF2 = 0;
countFailL2O = 0;
numRepts = 1;
for dmx = 1:numRepts
% Initialize Problem                                                        See testFiltersCappe.m for explaination of the problem initialization
problem = simple_benchmark();
initialState = 0;
measurementNoise = 2;
processNoise = 10;
%% Initialize UKF

ukf = unscentedKalmanFilter( ... 
@simple_ukfStateFunc,...
@simple_ukfDataGen, ...
initialState, ...
'HasAdditiveMeasurementNoise',1);
ukf.ProcessNoise = processNoise;
ukf.MeasurementNoise = measurementNoise;

ekf = extendedKalmanFilter( ...
@simple_ukfStateFunc,...
@simple_ukfDataGen, ...
initialState, ...
'HasAdditiveMeasurementNoise',1);
ekf.ProcessNoise = processNoise;
ekf.MeasurementNoise =measurementNoise;


%% Initiazile PF
numParticles = 1000;
% xx, yy, zz are x_t, x_t-1, time respectively
pdf_f = @(xx, yy, zz, dat) normpdf(xx, simple_benchmark.StateTransitionC(yy, zz), processNoise^0.5);
pdf_q = @(xx, yy, zz, dat) normpdf(xx, simple_benchmark.StateTransitionC(yy, zz), measurementNoise^0.5);
% xx, dat are current state and data respectively
pdf_g = @(dat,xx) normpdf(dat, simple_benchmark.GetDataC(xx), measurementNoise^0.5);

PF.resampling = 0;
PF.statePDF = pdf_f;
PF.likelihood = pdf_g;
PF.proposalPDF = pdf_q; % Bootstrap filter
%PF.proposalSampler.random = @(xx, dat, Tim) normrnd(simple_benchmark.StateTransitionC(xx, Tim), processNoise^0.5);
PF.proposalSampler.random = @(xx, dat, Tim) simple_benchmark.StateTransition(xx, Tim);
PF.particles = zeros(length(initialState), numParticles);
PF.weights = zeros(1, numParticles);
for idx = 1:numParticles
    data = simple_benchmark.GetDataC(initialState);
    PF.particles(:,idx) = PF.proposalSampler.random(initialState, ...
        data, 0);
    PF.weights(idx) = PF.likelihood(data, initialState);
end
PF.weights = PF.weights/sum(PF.weights);


%% L2O Filter
L2O.samples = linspace(-30, 500, 1000);
L2O.weights = ones(length(L2O.samples(1,:)), 1)./length(L2O.samples(1,:));
L2O.data = data;
L2O.likelihood = @(dat,outp) normpdf(dat, outp, 1);
L2O.stateEvolution = @(xx, yy) simple_benchmark.StateTransition(xx, yy); %xx is the state, yy is time;
L2O.oper = 0;
L2O.output = simple_benchmark.GetDataC(L2O.samples);
L2O.dataFunc = @(xx, yy) simple_benchmark.GetDataC(xx);
L2O.R = measurementNoise;
L2O.estimate = initialState;
L2O.err = processNoise;


%% Run Problem

state(1) = initialState;
    disp(dmx)
    clear err err2 err3 err4 err5
    clear est est2 est3 est4 est5
    clear state
    state(1) = initialState;
    for idx = 2:100
        % Evolve State and collect data
        state(idx) = simple_benchmark.StateTransition(state(idx - 1), idx);
        data = simple_benchmark.GetData(state(idx));
        % Evaluate UKF
        tic
        [est(idx), err(idx)] = correct(ukf, data);
        predict(ukf);
        UT(idx-1) = toc;

        % Evaluate Particle Filter
        PF.data = data;
        PF.oper = idx;
        tic
        PF = BootstrapFilter(PF);
        PT(idx-1) = toc;
        est2(idx) = PF.estimate;
        err2(idx) = PF.err;
        

        % Evaluate L2O_filter
        L2O.oper = idx;
        L2O.data = data;
        tic
        L2O = L2oFilter_up(L2O, 0, 1000, 500, 1);
        LT(idx-1) = toc;
        est3(idx) = L2O.estimate;
        err3(idx) = L2O.err;

        % Evaluate EKF
        tic
        [est4(idx), err4(idx)] = correct(ekf, data);
        predict(ekf);
        ET(idx-1) = toc;
        
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
    for tmz = 1:100
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

h = figure;
subplot(2, 2, 1);
errorbar(zeros(1,length(est)), 3*err.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est),length(est)), est-state, 'black', 'filled','LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 200])
legend('Error Bounds', 'True Error')
title('UKF')

subplot(2, 2, 2);
errorbar(zeros(1,length(est2)), 3*err2.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est2),length(est2)), est2-state, 'black', 'filled','LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 200])
title('Particle Filter')

subplot(2, 2, 3);
hold on
errorbar(zeros(1,length(est3)), 3*err3.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est3),length(est3)), est3-state, 'black','filled', 'LineWidth', 0.05)
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 200])
title('L2O Filter')

subplot(2, 2, 4);
errorbar(zeros(1,length(est4)), 3*err4.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est4),length(est4)), est4-state, 'black', 'filled', 'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('Error')
ylim([-60, 60])
xlim([0, 200])
title('EKF')


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

h = figure;
subplot(2, 2, 1);
plot(est, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
errorbar(est, 3*err.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est),length(est)), state, 'black', 'filled', 'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('State')
ylim([-20, 300])
xlim([0, 200])
legend('Estimate', 'Error Bounds', 'True State')
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
ylim([-20, 300])
xlim([0, 200])
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
ylim([-20, 300])
xlim([0, 200])
title('L2O Filter')
hold on

subplot(2, 2, 4);
plot(est4, 'blue', 'LineStyle', '-', 'LineWidth', 2.0)
hold on
errorbar(est4, 3*err4.^0.5, 'blue','LineStyle',  '-')
hold on
scatter(linspace(1,length(est4),length(est4)), state, 'black','filled',  'LineWidth', 0.05)
hold on
xlabel('Time')
ylabel('State')
ylim([-20, 300])
xlim([0, 200])
title('EKF')
hold on

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
hold on