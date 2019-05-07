
%% Case 1
figure(1)
% Uniform Load
X = 30*ones([100, 1001]);
X(:, end) = linspace(0, 1, 100);


[strains, M] = CantileverBeamTest(X);
disp(strains)

subplot(3, 1, 1)
plot(linspace(0,1, 1000), X(1, 1:end-1), 'k', 'LineWidth', 2.0)
hold on
ylabel('Load')
xlabel('X')
xlim([0, 1])
ylim([0, max(X(1, 1:end - 1)) + 5])

subplot(3, 1, 2)
plot(X(:, end), M, 'b', 'LineWidth', 2.0)
hold on
ylabel('Moment')
xlabel('X')
    
subplot(3, 1, 3)
plot(X(:, end), strains, 'b', 'LineWidth', 2.0)
hold on
ylabel('\epsilon')
xlabel('X')


%% Case 2
figure(2)
% Sinusoidal Load

X = 30*sin(linspace(0, 4*pi, 1001));
X = repmat(X, 100, 1);
X(:, end) = linspace(0, 1, 100);


[strains, M] = CantileverBeamTest(X);
disp(strains)

subplot(3, 1, 1)
plot(linspace(0,1, 1000), X(1, 1:end-1), 'k', 'LineWidth', 2.0)
hold on
ylabel('Load')
xlabel('X')
xlim([0, 1])
ylim([min(X(1, 1:end - 1)) - 5, max(X(1, 1:end - 1)) + 5])

subplot(3, 1, 2)
plot(X(:, end), M, 'b', 'LineWidth', 2.0)
hold on
ylabel('Moment')
xlabel('X')
    
subplot(3, 1, 3)
plot(X(:, end), strains, 'b', 'LineWidth', 2.0)
hold on
ylabel('\epsilon')
xlabel('X')


%% Case 2
figure(3)
% Triangle Load

X = 30*linspace(1, -1, 1001);
X = repmat(X, 100, 1);
X(:, end) = linspace(0, 1, 100);


[strains, M] = CantileverBeamTest(X);
disp(strains)

subplot(3, 1, 1)
plot(linspace(0,1, 1000), X(1, 1:end-1), 'k', 'LineWidth', 2.0)
hold on
ylabel('Load')
xlabel('X')
xlim([0, 1])
ylim([min(X(1, 1:end - 1)) - 5, max(X(1, 1:end - 1)) + 5])

subplot(3, 1, 2)
plot(X(:, end), M, 'b', 'LineWidth', 2.0)
hold on
ylabel('Moment')
xlabel('X')
    
subplot(3, 1, 3)
plot(X(:, end), strains, 'b', 'LineWidth', 2.0)
hold on
ylabel('\epsilon')
xlabel('X')
