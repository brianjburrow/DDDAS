clc
close all
clear all


nParams = 100;
paramSet = zeros([nParams, 2]);


a = 0.0005;
b = 0.01;
paramSet(:,1) = (b-a).*rand(nParams,1) + a;

a = 0.00005;
b = 0.0005;
paramSet(:,2) = (b-a).*rand(nParams,1) + a;

damage = zeros([nParams, 100, 2]);
for iParams = 1:nParams
    for iTime = 2:250
        damage(iParams, iTime, :) = evolve(damage(iParams, iTime - 1,:), paramSet(iParams,:));
    end
end
figure(1)
xlim([0,1])
hold on
ylim([0,1])
xlabel("Damage Parameter 1")
hold on
ylabel("Damage Parameter 2")
for iParams = 1:nParams
    scatter(damage(iParams, :,1), damage(iParams,:,2), 'filled', 'black');
    drawnow
    pause(0.1)
    hold on
end
function damage = evolve(damage, params)
    damage = damage(1,:,:);
    damage = reshape(damage, [1,2]);
    nCases = length(damage(:,1));
    for iCase = 1:nCases
        if damage(iCase, 1) < 1
            pd = makedist('Normal', 'mu', 0.0005, 'sigma', 0.005/(0.05*damage(iCase,1)+1));
            t = truncate(pd, 0, 1 - damage(iCase,1));
            damage(iCase, 1) = damage(iCase,1) + random(t, 1);
        end
        if damage(iCase, 2) < .99
            pd2 = makedist('Normal', 'mu', params(1,1)*damage(iCase, 1), 'sigma', params(1,2)*damage(iCase, 1) + 0.000001);
            t2 = truncate(pd2, 0, 1 - damage(iCase,2));
            damage(iCase, 2) = damage(iCase, 2) + random(t2, 1);
        end
    end
end