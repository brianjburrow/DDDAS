clear all
close all

polyOrders = [1, 2, 3, 4];
dim = 1;  % Dimensionality of parameter space

%samps = genPosteriorSamples(4);
samps = normrnd(100, 10, [200,1]);
K = length(samps(:,1));

options = optimoptions('fmincon', 'MaxIterations',  1000, 'algorithm',...
    'sqp', 'StepTolerance', 1e-14, 'MaxFunctionEvaluations', 10000);

for dmx = 1:length(polyOrders)
    for idx = 1:length(dim)
        multi_indices = genTotalOrderMI(polyOrders(dmx), idx);
        %multi_indices = genNoCross(polyOrders(dmx), idx);
        %multi_indices = genNoMixedMI(polyOrders(dmx), idx);
        disp(multi_indices)
        if dmx == 1
            gammas(dmx).val = ones([1, length(multi_indices(:,1))]);
        else
            disp(gammas(dmx-1).val)
            extras = length(multi_indices(:,1)) - length(gammas(dmx-1).val);
            gammas(dmx).val = [gammas(dmx-1).val, zeros([1, extras])];
            disp(gammas(dmx).val)
        end
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
    figure(dmx)
    r = zeros(size(samps));
    for idx = 1:length(samps(:,1))
        r(idx) = tMAP(samps(idx,:), gammas(dmx).val, multi_indices);
    end

    histogram(samps, 30)
    hold on
    histogram(r, 30, 'FaceColor', 'red')
end
