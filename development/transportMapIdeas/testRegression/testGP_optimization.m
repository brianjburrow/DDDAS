clear all
close all

polyOrder = 3;
dim = 1;  % Dimensionality of parameter space
lb = 1;
ub = 5;
numSlice = 10;
numSamples = 400;
%samps = genPosteriorSamples(4);


%% Run Script
xxx = linspace(lb, ub, numSlice);
samples = zeros(numSamples, numSlice);

for idx = 1:numSlice
    variance = 10*abs(sin(xxx(idx)));
    samples(:,idx) = 15 + betarnd(3, 2*xxx(idx), [numSamples,1]);
    scatter(xxx(idx)*ones([numSamples,1]), samples(:,idx), 'blue');
    hold on
    xlabel("Input")
    ylabel("P_{tar}")
end

K = length(samples(:,1));

options = optimoptions('fmincon', 'MaxIterations',  1000, 'algorithm',...
    'sqp', 'StepTolerance', eps, 'MaxFunctionEvaluations', 10000);

gammas_1 = zeros([numSlice, 1]);
gammas_2 = zeros([numSlice, 1]);
gammas_3 = zeros([numSlice, 1]);
gammas_4 = zeros([numSlice, 1]);

for tmz = 1:polyOrder
    fprintf("TMZ %d \n", tmz)
    multi_indices = genTotalOrderMI(tmz, 1);
    for dmx = 1:numSlice
        fprintf("DMX %d \n", dmx)
        samps = samples(:,dmx);
        if tmz == 1
            objF  = @(xx) objectiveFunc(xx, samps, multi_indices, 1);
            const = @(xx) constraintFunc(xx, samps, multi_indices, 1, 10^-8);
            if dmx == 1
                AA = fmincon(...
                    objF, ...
                    [1,1],...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
                gammas_1(dmx) = AA(1);
                gammas_2(dmx) = AA(2);
            else
                AA = fmincon(...
                    objF, ...
                    [gammas_1(dmx-1), gammas_2(dmx-1)],...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
                gammas_1(dmx) = AA(1);
                gammas_2(dmx) = AA(2);
            end
        elseif tmz == 2
            objF  = @(xx) objectiveFunc([gammas_1(dmx), gammas_2(dmx), xx], samps, multi_indices, 1);
            const = @(xx) constraintFunc([gammas_1(dmx),gammas_2(dmx),xx], samps, multi_indices, 1, 10^-8);
            if dmx == 1
                gammas_3(dmx) = fmincon(...
                    objF, ...
                    1,...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
            else
                gammas_3(dmx) = fmincon(...
                objF, ...
                gammas_3(dmx-1),...
                [],[],[],[],[],[],...
                const, ...
                options...
                );
            end
        elseif tmz == 3
            objF  = @(xx) objectiveFunc([gammas_1(dmx), gammas_2(dmx), gammas_3(dmx), xx], samps, multi_indices, 1);
            const = @(xx) constraintFunc([gammas_1(dmx), gammas_2(dmx), gammas_3(dmx), xx], samps, multi_indices, 1, 10^-8);
            if dmx == 1
                gammas_4(dmx) = fmincon(...
                    objF, ...
                    1,...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
            else
                gammas_4(dmx) = fmincon(...
                    objF, ...
                    gammas_4(dmx-1),...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
            end
        end
    end
end

transformedSamples = zeros(size(samples));

for idx = 1:numSlice
    for tmz = 1:numSamples
        transformedSamples(tmz,idx) = tMAP(samples(tmz,idx), [gammas_1(idx), gammas_2(idx), gammas_3(idx), gammas_4(idx)], multi_indices);
    end
    subplot(6,1,6)
    scatter(xxx(idx)*ones(size(samples(:,idx))), transformedSamples(:,idx), 'red') 
    ylabel("P_{ref}")
    hold on
end

map_regressor_1  = fitrgp(xxx', gammas_1  ,'KernelFunction','squaredexponential', 'Standardize',true, 'OptimizeHyperparameters','auto');
map_regressor_2  = fitrgp(xxx', gammas_2  ,'KernelFunction','squaredexponential', 'Standardize',true, 'OptimizeHyperparameters','auto');
map_regressor_3  = fitrgp(xxx', gammas_3  ,'KernelFunction','squaredexponential', 'Standardize',true, 'OptimizeHyperparameters','auto');
map_regressor_4  = fitrgp(xxx', gammas_4  ,'KernelFunction','squaredexponential', 'Standardize',true, 'OptimizeHyperparameters','auto');

close all
xms = linspace(lb, ub, 100)';
[y1, e1] = predict(map_regressor_1, xms);
[y2, e2] = predict(map_regressor_2, xms);
[y3, e3] = predict(map_regressor_3, xms);
[y4, e4] = predict(map_regressor_4, xms);

plot3dhist(samples, numSlice, xxx, 13, false)

plot3dhist(transformedSamples, numSlice, xxx, 14, true)

figure(100)
subplot(6,1,1)
for idx = 1:numSlice
    scatter(xxx(idx)*ones([numSamples,1]), samples(:,idx), 'blue');
    hold on
    xlabel("Input")
    ylabel("P_{tar}")
end
set(gca,'xtick',[])

figure(1000)
subplot(4,1,1)
scatter(xxx, gammas_1)
hold on
errorbar(xms, y1, 3*e1)
ylabel('P 1')
set(gca,'xtick',[])
subplot(4,1,2)
scatter(xxx, gammas_2)
hold on
errorbar(xms, y2, 3*e2)
ylabel('P 2')
set(gca,'xtick',[])
subplot(4,1,3)
scatter(xxx, gammas_3)
hold on
errorbar(xms, y3, 3*e3)
ylabel('P 3')
set(gca,'xtick',[])
subplot(4,1,4)
scatter(xxx, gammas_4)
hold on
errorbar(xms, y4, 3*e4)
ylabel('P 3')



figure(1203431)
subplot(6,1,6)
for idx = 1:numSlice
    scatter(xxx(idx) * ones([numSamples,1]), transformedSamples(:,idx), 'red')
    hold on
    xlabel('Input')
    ylabel("P_{ref}")
end


