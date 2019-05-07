% Finding optimal importance weights - S Amaral, D Allaire, K Willcox
% Rewritten in MATLAB - B Isaac

clear all
clc
format long
rng(1)
% Initialize
n=100;      %// number of proposal samples
m=1000;       %// number of target samples              
d=1;          %// number of dimensions
itr=1000;     %// number of frank-wolfe iterations
X = zeros(n,d);
Y = zeros(m,d);
w = 1/n .* ones(n,1);
feval = zeros(itr+1,1);

% Proposal distribution (Multivariate Gaussian)
muX = zeros(d,1);
covX = zeros(d,d);
for i=1:d
    muX(i) = 1/sqrt(d);
    covX(i,i) = 2;
end

% covXL = zeros(d,d);
covXL = chol(covX);

%// Target distribution (Multivariate Gaussian)
muY = zeros(d,1);
covY = zeros(d,d);
for i=1:d
    covY(i,i) = 0.025;
end

% covYL = zeros(d,d);
covYL = chol(covY);

%  Generate (nXd) independent Gaussian RV's
%t = normrnd(0,1,n,1);
t = 0;

% Proposal samples
for i=1:n
    if i > 0
        t = t + 1/n;
        disp(t)
    end
    for j=1:d
          for k=1:d
                X(i,j) = X(i,j)+ covXL(j,k)*t(k) + muX(j); 
          end
    end
end

X(:, 1) = linspace(-100,-20,n);

% Target samples
for i=1:m
    t = normrnd(0,1,d,1);
    for j=1:d
        for k=1:d
           Y(i,j) = Y(i,j) + covYL(j,k)*t(k) + muY(j); 
        end
    end
end
Y = Y.^3;
%  Resize into unit-hypercube
[X,Y] = unitHyper(X,Y);

% Compute matrix-vector multiply (H * w)
Hw = matvecProd(X,w);

%  Compute the RHS vector b
b = buildXYb(X,Y);

%  Solve the optimization statement using the Frank-Wolfe algorithm
[w, feval] = frankwolfe(X,b, w, Hw, feval, itr);
disp('done')
d=1;
% Plot results EDF
[M1,idx1] = sort(X(:,d));
[M2,idx2] = sort(Y(:,d));
 
figure(1)
% Proposal
stairs(M1,(1/n:1/n:1),'r', 'LineStyle', '-', 'LineWidth', 2.0)
xlabel('Value')
ylabel('Distribution Function');
ylim([0,1])
legend('Proposal ECDF')
hold on
% Target
stairs(M2,(1/m:1/m:1),'b', 'LineStyle', ':', 'LineWidth', 2.0)
legend('Proposal ECDF','Target ECDF')
hold on
% Weighted Proposal
wnew = zeros(n,1);
wnew(1)=0;
for i=1:n
    wnew(i+1) = wnew(i)+w(idx1(i));
end

stairs([0;M1],wnew,'k', 'LineStyle', '--', 'LineWidth', 2.0);
xlabel('Value')
ylabel('Distribution Function');
ylim([0,1])
legend('Proposal ECDF','Target ECDF', 'Weighted Proposal ECDF'); 

% Plot results Convergence
%void plotConv(feval, m2, duration,option, d) 
option = 1;
m2 = Y;
m   = size(m2,1);
itr = size(feval,1);
minC = min(feval);
maxC = max(feval);

y = -minC+0.001;
if (option == 1) 
    %v = 1/m*ones(m,1);
    v = buildXYb(m2,m2);
    for i=1:m
      y = y+ 0.5*v(i)/m;
    end
end 

figure(2)
for ii=1:itr
    xf(ii) = (ii-1)/(itr);
end
plot(xf,feval+y);
xlabel('Computational Time [sec]');
ylabel('log_{10end(Objective Function)}');
legend('Proposal'); 



