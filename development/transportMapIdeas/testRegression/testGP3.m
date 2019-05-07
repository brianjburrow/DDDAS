clear all
close all

polyOrders = 2;
dim = 1;  % Dimensionality of parameter space
lb = 1;
ub = 5;
numSlice = 50;
numSamples = 400;
%samps = genPosteriorSamples(4);


%% Run Script
xxx = linspace(lb, ub, numSlice);
samples = zeros(numSamples, numSlice);
figure(1)
subplot(5,1,1)
for idx = 1:numSlice
    variance = 10*abs(sin(xxx(idx)));
    samples(:,idx) = betarnd(xxx(idx), 5, [numSamples,1]);
    scatter(xxx(idx)*ones([numSamples,1]), samples(:,idx), 'blue');
    hold on
    xlabel("Input")
    ylabel("P_{tar}")
end


samps = normrnd(100, 10, [200,1]);
K = length(samps(:,1));

options = optimoptions('fmincon', 'MaxIterations',  1000, 'algorithm',...
    'sqp', 'StepTolerance', eps, 'MaxFunctionEvaluations', 10000);


figure(1)
subplot(5,1,2)
xlim([1, 10])
ylim([0, 100])
plot([0,5], [2,2],'blue')
hold on
plot([0,5], [-2,-2], 'blue')
hold on
for dmx = 1:numSlice
    fprintf("Iteration %d", dmx);
    samps = samples(:,dmx);
    for idx = 1:length(dim)
        multi_indices = genTotalOrderMI(polyOrders, idx);
        %multi_indices = genNoCross(polyOrders(dmx), idx);
        %multi_indices = genNoMixedMI(polyOrders(dmx), idx);
        if dmx == 1
            gammas(dmx).val = ones([1, length(multi_indices(:,1))]);
        else
            gammas(dmx).val = gammas(dmx-1).val;
        end
        if dmx == 1
            for kbtx = 1:100
                disp("Proposal")
                randInit = unifrnd(0, 1000, size(gammas(dmx).val));
                disp(randInit)
                objF  = @(xx) objectiveFunc(xx, samps, multi_indices, idx);
                const = @(xx) constraintFunc(xx, samps, multi_indices, idx, 10^-8);
                [tempOpt, tempObj] = fmincon(...
                            objF, ...
                            randInit,...
                            [],[],[],[],[],[],...
                            const, ...
                            options...
                            );
                if kbtx > 1 
                    if tempObj <= prevObj
                        gammas(dmx).val = tempOpt;
                    end
                else
                    gammas(dmx).val = tempOpt;
                end
                disp("Found")
                disp(tempOpt)
                disp("Selected")
                disp(gammas(dmx).val)
                prevObj = tempObj;
            end
        else
            objF  = @(xx) objectiveFunc(xx, samps, multi_indices, idx);
            const = @(xx) constraintFunc(xx, samps, multi_indices, idx, 10^-8);
            gammas(dmx).val = fmincon(...
                objF, ...
                gammas(dmx).val,...
                [],[],[],[],[],[],...
                const, ...
                options...
                );
        end
    end
    r = zeros(size(samps));
    for idx = 1:length(samps(:,1))
        r(idx) = tMAP(samps(idx,:), gammas(dmx).val, multi_indices);
    end

    figure(1)
    subplot(5,1,2)
    xlim([1, 5])
    plot([0,5], [2,2],'blue')
    hold on
    plot([0,5], [-2,-2], 'blue')
    hold on
    scatter(xxx(dmx)*ones([numSamples,1]), r, 'red');
    hold on
    xlabel("Input")
    ylabel("P_{ref}")
    drawnow 
    
    figure(1)
    subplot(5,1,3)
    xlim([1, 5])
    hold on
    scatter(xxx(dmx), gammas(dmx).val(1), 'black');
    hold on
    xlabel("Input")
    ylabel({"MAP", "Param 1"})
    drawnow 
    
    figure(1)
    subplot(5,1,4)
    xlim([1, 5])
    hold on
    scatter(xxx(dmx), gammas(dmx).val(2), 'black');
    hold on
    xlabel("Input")
    ylabel({"MAP", " Param 2"})
    drawnow 
    
    figure(1)
    subplot(5,1,5)
    xlim([1, 5])
    hold on
    scatter(xxx(dmx), gammas(dmx).val(3), 'black');
    hold on
    xlabel("Input")
    ylabel({"MAP", " Param 3"})
    drawnow 
end

for tmz = 1:numSlice
    x1(tmz) = gammas(tmz).val(1);
    x2(tmz) = gammas(tmz).val(2);
    x3(tmz) = gammas(tmz).val(3);
end

x1r = fitrsvm(xxx', x1,'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true, 'OptimizeHyperparameters','auto');

x2r = fitrsvm(xxx', x2,'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true, 'OptimizeHyperparameters','auto');

x3r = fitrsvm(xxx', x3,'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true, 'OptimizeHyperparameters','auto');


figure(1)
subplot(5,1,3)
plot(xxx, predict(x1r, xxx'), 'red')
drawnow 

figure(1)
subplot(5,1,4)
plot(xxx, predict(x2r, xxx'), 'red')


figure(1)
subplot(5,1,5)
plot(xxx, predict(x3r, xxx'), 'red')