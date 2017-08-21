clear

%% System model
% Plant matrices
A = [0,1;0,0];
b = [0;1];
d = length(b); %system size
% initial states
x0 = [1;1];
% Horizon length
T = 5;

%% Time discretization
% Discretization size
n = 10000; % grid size
h = T/n; % discretization interval
% System discretization
[Ad,bd] = c2d(A,b,h);
% Matrix Phi
Phi = zeros(d,n);
v = bd;
Phi(:,end) = v;
for j = 1:n-1
    v = Ad*v;
    Phi(:,end-j) = v;
end
% Vector zeta
zeta = -Ad^n*x0;

%% Convex optimizaiton via CVX
cvx_begin
 variable u(n)
 minimize norm(u,1)
 subject to 
   Phi*u == zeta;
   norm(u,inf) <= 1;
cvx_end

%% Plot
figure;
plot(0:T/n:T-T/n,u);
title('Sparse control');


