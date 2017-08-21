
clear 
%% polynomial coefficients
x_orig = [-1,zeros(1,78),1,0]';
%% sampling
t = 0:0.1:1;
y = polyval(x_orig,t);
%% data size
N = length(t);
M = N-1;
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
%% Plot
tt = 0:0.01:1;
figure;
stem(t,y); hold on
plot(tt,polyval(x_orig,tt),'--');
plot(tt,polyval(x,tt));
