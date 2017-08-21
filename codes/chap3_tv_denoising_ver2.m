clear;

%% Read image
X_orig=imread('lena.jpg');
%X_orig = phantom('Modified Shepp-Logan',200);
[n,m]=size(X_orig);

%% Noise
rng(1);
Y = X_orig + uint8(randi(50,n,m)); % for lena (1)
%Y = X_orig + uint8(255*(randi(2,n,m)-1)); %for lena (2)
%Y = X_orig + rand(n,m)*0.5; % for phantom

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
 Z = zeros(n,m); V = zeros(n,m); % initial values
 M = Phi'*Phi + (1/gamma)*Psi'*Psi; % a matrix in the first step

Yv = double(Y);
W = Phi'*Yv;
for k=1:N
    X = M\(W+gamma\Psi'*(Z-V));
        P = Psi*X+V;
        Z = soft_thresholding(gamma*lambda,P);
        V = P - Z;
end

%% Result
figure;
imshow(uint8(round(X))); % for lena
%imshow(X); % for phantom
title('Restored image');




