%% Case 2
figure(2)
% Sinusoidal Load

X = [30*sin(linspace(0, 1.5*pi, 500)) + 100*ones([1, 500]), 15*sin(linspace(1.5*pi/4, 1.75*pi, 501))+ 55*ones([1, 501])];
X = repmat(X, 100, 1);
X(:, end) = linspace(0, 1, 100);

X = computeNoisyLoads(X);

[strains, M] = CantileverBeamTest(X);

subplot(3, 1, 1)
for iSample = 1:100
    plot(linspace(0, 1, 1000), X(iSample, 1:1000), 'k', 'LineWidth', 2.0)
    hold on
    xlabel('X')
    ylabel('Loads')
    xlim([0, 1])
    ylim([0, 150])
    hold off
    pause(0.005)
end

subplot(3, 1, 2)
histogram(M)
hold on
xlabel("X")
ylabel("Moment")

subplot(3, 1, 3)
histogram(strains)
hold on
xlabel("X")
ylabel("Strains")

function loads = computeNoisyLoads(X)
    [nLoads, nElem] = size(X);
    loads = zeros([nLoads, nElem]);
    for iLoad = 1:nLoads
        loads(iLoad, :) = [...
            X(iLoad, 1:end-1) + mvnrnd(zeros([1, 1000]), 0.5*eye(1000)), ...
        X(iLoad, end)];
    end
end
