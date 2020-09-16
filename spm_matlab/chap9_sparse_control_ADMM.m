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
n = 1000; % grid size
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

%% Convex optimizaiton via ADMM
mu = 2*n+d;
Psi = [eye(n);eye(n);Phi];
%M = (Psi'*Psi)\Psi';
M = (0.5*eye(n) - 0.5*Phi'*inv(2*eye(d)+Phi*Phi')*Phi)*Psi';
sat = @(x) sign(x).*min(abs(x),1);

EPS = 1e-5; % if the residue < EPS then the iteration will stop
MAX_ITER = 100000; % maximum number of iterations
z = [zeros(2*n,1);zeta]; v = zeros(mu,1);
r = zeta;
k = 0;
gamma = 0.05;

tic
while (norm(r)>EPS) & (k < MAX_ITER)
    u = M*(z-v);
    z0 = soft_thresholding(gamma,u+v(1:n));
    z1 = sat(u+v(n+1:2*n));
    z2 = zeta;
    z = [z0;z1;z2];
    v = v + Psi*u - z;
    r = Phi*u - zeta;
    k = k + 1;
end
cpt=toc;

%% Plot
figure;
plot(0:T/n:T-T/n,u,'LineWidth',2);
title('Sparse control');


