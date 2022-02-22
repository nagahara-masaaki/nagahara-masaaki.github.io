%% An example of over fitting
%% for IEICE review article
clear
rng(1)

%% data
t = 0:0.1:1.5;
y_orig = 2*t;

%% noise
ns = randn([length(t),1])*0.1; % noise (mean=0, std deviation=0.1)
y = y_orig(:)+ns(:);

%% interpolation
A = vander(t);
x = inv(A)*y;

%% regularization
lambda = 1;
xr = inv(lambda*eye(length(t))+A'*A)*A'*y;


%% LASSO by FISTA
lambda_s = 0.5;

%% FISTA algorithm
n = length(t);
N = 100000; % max number of iterations
z0 = xr; % initial guess for z
x0 = z0; % initial guess for x
tau0 = 1; % initial guess for tau

gamma = 1/norm(A,2)^2; % parameter of proximal operator

X1 = x0;
Z1 = z0;
tau1 = tau0;
for k = 1:N
    X2 = soft_thresholding(gamma*lambda_s,Z1-gamma*A'*(A*Z1-y));
    tau2 = (1+sqrt(1+4*tau1^2))/2;
    Z2 = X2 + (tau1-1)/tau2*(X2 - X1);
    X1 = X2;
    Z1 = Z2;
    tau1 = tau2;
end

xs = Z1;

%% L0 regularization by IHT
lambda_s = 0.5;

%% IHT algorithm
n = length(t);
N = 1000000; % max number of iterations
x0 = xs; % initial guess

gamma = 1/norm(A,2)^2; % parameter of proximal operator

X = x0;
for k=1:N
    X = hard_thresholding(sqrt(2*gamma*lambda_s),X-gamma*A'*(A*X-y));
end

xs_polished = X;

%% plot
tt = t(1):0.01:t(end);
yy = polyval(x,tt);
yyr = polyval(xr,tt);
yys = polyval(xs,tt);

figure;
plot(tt,yy,'LineWidth',1);
hold on
stem(t,y,'LineWidth',1,'MarkerSize',10,'MarkerFaceColor','y');
grid on
xlabel('t')
ylabel('y')
set(gca,'FontSize',18)

figure;plot(tt,yyr,'LineWidth',1);
hold on
stem(t,y,'LineWidth',1,'MarkerSize',10,'MarkerFaceColor','y');
%plot(tt,yy,'b--','LineWidth',0.5);
grid on
xlabel('t')
ylabel('y')
set(gca,'FontSize',18)

figure;plot(tt,yys,'LineWidth',1);
hold on
stem(t,y,'LineWidth',1,'MarkerSize',10,'MarkerFaceColor','y');
plot(tt,yyr,'b--','LineWidth',0.5);
grid on
xlabel('t')
ylabel('y')
set(gca,'FontSize',18)

disp('LASSO solution')
xs

disp('Polished solution by IHT')
xs_polished

function sv = soft_thresholding(lambda,v)

[m,n]=size(v);
mn = m*n;
sv = zeros(m,n);
for i = 1:mn
    if abs(v(i))<=lambda
        sv(i) = 0;
    else
        sv(i) = v(i) - sign(v(i))*lambda;
    end
end

end

function hv = hard_thresholding(lambda,v)

[m,n]=size(v);
mn = m*n;
hv = zeros(m,n);
for i = 1:mn
    if abs(v(i))<=lambda
        hv(i) = 0;
    else
        hv(i) = v(i);
    end
end

end
