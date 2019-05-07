clear all
close all

polyOrder = 4;
dim = 1;  % Dimensionality of parameter space
lb = 1;
ub = 5;
numSlice = 1;
numSamples = 20000;
burnIn = 10000;

%% Run Script
xxx = linspace(lb, ub, numSlice);
samples = zeros(numSamples, numSlice);
figure(1)
subplot(7,1,1)
for idx = 1:numSlice
    variance = 10*abs(sin(xxx(idx)));
    %samples(:,idx) = genPosteriorSamples(xxx(idx), numSamples, burnIn);
    scatter(xxx(idx)*ones([numSamples,1]), samples(:,idx), 'blue');
    hold on
    xlabel("Input")
    ylabel("P_{tar}")
end


K = length(samples(:,1));



gammas_1 = zeros([numSlice, 1]);
gammas_2 = zeros([numSlice, 1]);
gammas_3 = zeros([numSlice, 1]);
gammas_4 = zeros([numSlice, 1]);
gammas_5 = zeros([numSlice, 1]);

for tmz = 1:polyOrder
    fprintf("TMZ %d \n", tmz)
    multi_indices = genTotalOrderMI(tmz, 1);
    
    if tmz > 1
        options = optimoptions('fmincon', 'MaxIterations',  1000, 'algorithm',...
                'sqp', 'StepTolerance', eps, 'MaxFunctionEvaluations', 10000,'UseParallel',true);
    else
        options = optimoptions('fmincon', 'MaxIterations',  1000, 'algorithm',...
                'sqp', 'StepTolerance', eps, 'MaxFunctionEvaluations', 10000,'UseParallel',false);
    end
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
        elseif tmz == 4
            objF  = @(xx) objectiveFunc([gammas_1(dmx), gammas_2(dmx), gammas_3(dmx), gammas_4(dmx), xx], samps, multi_indices, 1);
            const = @(xx) constraintFunc([gammas_1(dmx), gammas_2(dmx), gammas_3(dmx), gammas_4(dmx), xx], samps, multi_indices, 1, 10^-8);
            if dmx == 1
                gammas_5(dmx) = fmincon(...
                    objF, ...
                    1,...
                    [],[],[],[],[],[],...
                    const, ...
                    options...
                    );
            else
                gammas_5(dmx) = fmincon(...
                    objF, ...
                    gammas_5(dmx-1),...
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
        transformedSamples(tmz,idx) = tMAP(samples(tmz,idx), [gammas_1(idx), gammas_2(idx), gammas_3(idx), gammas_4(idx), gammas_5(idx)], multi_indices);
    end
    subplot(7,1,7)
    scatter(xxx(idx)*ones(size(samples(:,idx))), transformedSamples(:,idx), 'red') 
    ylabel("P_{ref}")
    hold on
end

subplot(7,1,2)
scatter(1:numSlice, gammas_1)
ylabel('P 1')
subplot(7,1,3)
scatter(1:numSlice, gammas_2)
ylabel('P 2')
subplot(7,1,4)
scatter(1:numSlice, gammas_3)
ylabel('P 3')
subplot(7,1,5)
scatter(1:numSlice, gammas_4)
ylabel('P 4')
subplot(7,1,6)
scatter(1:numSlice, gammas_4)
ylabel('P 5')



plot3dhist(samples, numSlice, xxx, 13, false)

plot3dhist(transformedSamples, numSlice, xxx, 14, true)