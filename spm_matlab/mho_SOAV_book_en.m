% SOAV control 

clear;
s = tf('s');


% Horizon length
T = 10;
N = 1000;
h = T/N;

% Plant
Ac = [0,1,0,0;0,0,1,0;0,0,0,1;0,0,0,0];
Bc = [0;0;0;1];
Cc = eye(4);
Dc = 0;

[n,m] = size(Bc);
Pss = ss(Ac,Bc,eye(n,n),zeros(n,m));
Pds = c2d(Pss,h);
[A,B,C,D] = ssdata(Pds);
n = size(A,1);

% Initial and terminal states
x0 = ones(n,1)/2;
xf = zeros(n,1);

% Discrete-valued control
Umin = -1;
Umax = 1;

alphabet = Umin:0.5:Umax;
p = ones(1,length(alphabet));
Np = length(p);

% L(u)
a = zeros(Np+1,1);
b = zeros(Np+1,1);
rm = alphabet*p';
a(1) = -1; a(Np+1) = 1;
b(1) = rm; b(Np+1) = -rm;

for k=2:Np
    a(k)=sum(p(1:k-1))-sum(p(k:Np));
    b(k)=-sum(p(1:k-1).*alphabet(1:k-1))+sum(p(k:Np).*alphabet(k:Np));
end

figure;
ug = Umin:0.01:Umax;
for i=1:Np
    plot(ug,p(i)*abs(ug-alphabet(i)),'--');
    hold on;
end
plot(ug,max(kron(ug,a)+kron(ones(1,length(ug)),b)),'r');
title('L(u)');
xlabel('u');
ylabel('L(u)');

% System equation
[Phi,Ups] = get_mat2(A,B,N);
V = [zeros(n,(N-1)*n),eye(n,n)];
g = V*Ups*x0;
F = V*Phi;

% piecewise linear cost function

AA = kron(eye(N),a);
BB = kron(ones(N,1),b);
EE = kron(eye(N),ones(length(a),1));

% Optimization
cvx_begin
variable U(N)
variable theta(N)
minimize(sum(theta))
subject to
g+F*U==xf;
U>=Umin;
U<=Umax;
%norm(U,inf)<=Umax;
%U>=0; % for non-negative control
AA*U+BB<=EE*theta;
cvx_end


% Optimization
cvx_begin
variable U1(N)
minimize(norm(U1,1))
subject to
g+F*U1==xf;
norm(U1,inf)<=Umax;
cvx_end

% Optimization
cvx_begin
variable U2(N)
minimize( norm(U2,2) )
subject to
g+F*U2==xf;
norm(U2,inf)<=Umax;
cvx_end

figure;
plot(0:h:h*(N-1),U,'r');
hold on
plot(0:h:h*(N-1),U1,'k-.');

plot(0:h:h*(N-1),U2,'--');
xlabel('time (sec)');
ylabel('u(t)');
title('Optimal Control');
legend('SOAV','L^1 optimal','L^2 optimal');

if n==2
    X = Ups*x0 + Phi*U;
    X_1 = X(2:2:n*N);
    X_2 = X(1:2:n*N);
    
    X1 = Ups*x0 + Phi*U1;
    X1_1 = X1(2:2:n*N);
    X1_2 = X1(1:2:n*N);

    X2 = Ups*x0 + Phi*U2;
    X2_1 = X2(2:2:n*N);
    X2_2 = X2(1:2:n*N);
    
    figure;
    plot(X_1,X_2,'r');
    hold on
    plot(X1_1,X1_2,'k-.');
    plot(X2_1,X2_2,'--');
    xlabel('x1')
    ylabel('x2')
    title('state-space trajectory');
    legend('SOAV','L^1 optimal','L^2 optimal');
else
    Xe = zeros(N,n);
    X2e = zeros(N,n);
    
    X = Ups*x0 + Phi*U;
    X2 = Ups*x0 + Phi*U2;
    
    for j = 1:n
        Xe(:,j) = X(j:n:end)
        X2e(:,j) = X2(j:n:end)
    end
    
    figure;
    plot(0:h:h*(N-1),Xe);
    hold on
    plot(0:h:h*(N-1),U,'--');
    xlabel('time (sec)')
    ylabel('x_i(t)')
    title('state variables x_i(t) and control u(t)');
end

return;

length(U(U~=0))*h
