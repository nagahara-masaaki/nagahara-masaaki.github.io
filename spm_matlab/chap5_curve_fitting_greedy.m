clear;

%% data
% polynomial coefficients
x_orig = [-1,zeros(1,78),1,0]';

% sampling
t = 0:0.1:1;
y = polyval(x_orig,t)';

% data size
N = length(t);
M = N-1;

% Order of polynomial
M_l = length(x_orig)-1;

% Vandermonde matrix
Phi=[];
for m=0:M_l
 Phi = [t'.^m,Phi];
end

%% Sparse modeling
% iteration parameters
EPS=1e-5; % if the residue < EPS then the iteration will stop
MAX_ITER=100000; % maximum number of iterations

% L1 by CVX
cvx_begin
 variable x_l1(M_l+1)
 minimize norm(x_l1,1)
 subject to 
   Phi*x_l1 == y
cvx_end

% Matching Pursuit
[x_mp,nitr_mp]=MP(y,Phi,EPS,MAX_ITER);

% OMP
[x_omp,nitr_omp]=OMP(y,Phi,EPS,MAX_ITER);

% CoSaMP
s = length(supp(x_orig));
[x_cosamp,nitr_cosamp]=CoSaMP(y,Phi,s,EPS,MAX_ITER);

% IHT
lambda=0.001;
gamma=0.01;
[x_iht,nitr_iht]=IHT(y,Phi,lambda,gamma,EPS,MAX_ITER);

% iterative s-sparse
gamma=0.01;
[x_iss,nitr_iss]=iterative_s_sparse(y,Phi,s,gamma,EPS,MAX_ITER);

%% Analysis
% coefficient plot
figure;
subplot(2,3,1)
stem(1:M_l+1,x_l1,'filled');
axis([0,M_l+1,-1,1])
title('L1-OPT')

subplot(2,3,2)
stem(1:M_l+1,x_mp,'filled');
axis([0,M_l+1,-1,1])
title('MP')

subplot(2,3,3)
stem(1:M_l+1,x_omp,'filled');
axis([0,M_l+1,-1,1])
title('OMP')

subplot(2,3,4)
stem(1:M_l+1,x_iht,'filled');
axis([0,M_l+1,-1,1])
title('IHT')

subplot(2,3,5)
stem(1:M_l+1,x_iss,'filled');
axis([0,M_l+1,-1,1])
title('ISS')

subplot(2,3,6)
stem(1:M_l+1,x_cosamp,'filled');
axis([0,M_l+1,-1,1])
title('COSAMP')

% residue
r_l1 = y-Phi*x_l1; %residue
e_l1 = norm(r_l1,2); %l2 norm of residue

r_mp = y-Phi*x_mp; %residue
e_mp = norm(r_mp,2); %l2 norm of residue

r_omp = y-Phi*x_omp; %residue
e_omp = norm(r_omp,2); %l2 norm of residue

r_iht = y-Phi*x_iht; %residue
e_iht = norm(r_iht,2); %l2 norm of residue

r_iss = y-Phi*x_iss; %residue
e_iss = norm(r_iss,2); %l2 norm of residue

r_cosamp = y-Phi*x_cosamp; %residue
e_cosamp = norm(r_cosamp,2); %l2 norm of residue

error = [
    e_l1;
    e_mp;
    e_omp;
    e_iht;
    e_iss;
    e_cosamp
    ];

n_itr = [
    10;
    nitr_mp;
    nitr_omp;
    nitr_iht;
    nitr_iss;
    nitr_cosamp
    ];

return;
 
%% Plot
tt = 0:0.01:1;
% Interpolation polynomial with order 10
figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x_i,tt));
ylim([0,1]);
title('interpolating polynomial (ord=10)');


figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x_2r,tt));
ylim([0,1]);
title('L2 regularization (ord=80)');

figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x,tt));
ylim([0,1]);
title('sparse modeling (ord=80)');


return;
