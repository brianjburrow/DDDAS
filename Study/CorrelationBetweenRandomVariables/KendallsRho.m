%% Generate Correlated Gaussian samples

x = mvnrnd([0, 0], [1., 0.5; 0.5, 1.], 100000);

[f, xi] = ksdensity(x);

[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);

scatter(x(:,1), x(:,2), 1, 'b', 'filled')
hold on
[C,h] = contour(X,Y,Z, 20, 'LineWidth', 2.0);
colormap(hot)
xlabel('X1')
ylabel('X2')


%% Convert to uniform marginals through the standard mvn CDF function
x_unif = normcdf(x);

figure(2)
scatterhist(x_unif(:,1), x_unif(:,2), 'MarkerSize', 0.5, 'Color', 'b')
ax=get(gca,'children');
set(ax,'Marker','.')
xlabel('Y1')
ylabel('Y2')


%% Use Gaussian Copula to create correlated Beta samples
x_gumbel = evinv(x_unif(:,1),0,1.);
x_beta = betainv(x_unif(:,2), 10, 2);

figure(3)
scatterhist(x_gumbel, x_beta, 'MarkerSize', 0.5, 'Color', 'b')
ax=get(gca,'children');
set(ax,'Marker','.')
xlabel('Y1')
ylabel('Y2')

%% Fit Kernel Density
figure(4)
x = [x_gumbel, x_beta];
[f, xi] = ksdensity(x);

[X, Y] = meshgrid(xi(:,1), xi(:,2));
Z = griddata(xi(:,1), xi(:,2), f, X, Y);

scatter(x(:,1), x(:,2), 1, 'b', 'filled')
hold on
[C,h] = contour(X,Y,Z, 20, 'LineWidth', 2.0);
colormap(hot)
xlabel('Maximum River Level')
ylabel('Probability of Flooding')

[rho, pval] = corr(x, 'Type', 'Kendall');
