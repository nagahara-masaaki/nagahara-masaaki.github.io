clear 

%% polynomial coefficients
x_orig = [-1,zeros(1,78),1,0]';
%% sampling
t = 0:0.1:1;
y = polyval(x_orig,t);
%% data size
%N = length(t);
%M = N-1;
%% Vandermonde matrix
Phi_v = vander(t);
%% Interpolation polynomial with order 10
x_i = inv(Phi_v)*y';

%% LASSO
% Order of polynomial
M_l = 80;
% Design matrix
Phi_l=[];
for m=0:M_l
 Phi_l = [t'.^m,Phi_l];
end
% CVX
cvx_begin
 variable x(M_l+1)
 minimize norm(x,1)
 subject to 
   Phi_l*x == y'
cvx_end

%% L2 regularization
r = 0.1;
N = length(x_orig);
x_2r = inv(r*eye(N) + Phi_l'*Phi_l)*Phi_l'*y';


%% Plot
tt = 0:0.01:1;
% Interpolation polynomial with order 10
figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x_i,tt),'LineWidth',2);
ylim([0,1]);
title('interpolating polynomial (ord=10)');
figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x_2r,tt),'LineWidth',2);
ylim([0,1]);
title('L2 regularization (ord=80)');
figure;
stem(t,y);
hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x,tt),'LineWidth',2);
ylim([0,1]);
title('sparse modeling (ord=80)');
