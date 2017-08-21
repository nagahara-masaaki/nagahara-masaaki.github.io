clear;

%% Read image
X_orig=imread('lena.jpg');
[n,m]=size(X_orig);

%% Noise
rng(1);
Y = X_orig + uint8(randi(50,n,m));

%% Display images
figure;
imshow(X_orig);
title('Original image');
figure;
imshow(Y);
title('Noisy image');

%% Denoising
% optimization parameter
 lambda = 50;
 
 Phi = eye(n);
 Psi = -diag(ones(n,1))+diag(ones(n-1,1),1);
 
 % ADMM iteration
 gamma = 1; % step size parameter
 N = 100; % number of iterations
 X_res = zeros(n,m); % restored image
 z = zeros(n,1); v = zeros(n,1); % initial values
 M = Phi'*Phi + (1/gamma)*Psi'*Psi; % a matrix in the first step

 for i=1:m
    %y = double(X_orig(:,i));
    y = double(Y(:,i));
    w = Phi'*y;
    for k=1:N
        x = M\(w+gamma\Psi'*(z-v));
        p = Psi*x+v;
        z = soft_thresholding(gamma*lambda,p);
        v = p - z;
    end
    X_res(:,i)=x;
end

%% Result
figure;
imshow(uint8(round(X_res)));
title('Restored image');




