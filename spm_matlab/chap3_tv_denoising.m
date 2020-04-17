clear;

%% Read image
Img = imread('cat.jpg'); % read image
X_orig = rgb2gray(Img); % Color to gray
[n,m] = size(X_orig); % Image size

%% Noise (salt & pepper)
rng(1);
Y = imnoise(X_orig,'salt & pepper',0.05);

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

% matrix Phi and Psi
Phi = eye(n);
Psi = -diag(ones(n,1))+diag(ones(n-1,1),1);

% ADMM iteration
gamma = 1; % step size parameter
N = 500; % number of iterations
X_res = zeros(n,m); % restored image
Z = zeros(n,m); V = zeros(n,m); % initial values
M = sparse(Phi'*Phi + (1/gamma)*Psi'*Psi); % matrix in the first step (sparse matrix)

% vertical processing
disp('optimization may take a lot of time...')
Yv = double(Y);
W = Phi'*Yv;
for k=1:N
    X = M\(W+gamma\Psi'*(Z-V));
    P = Psi*X+V;
    Z = soft_thresholding(gamma*lambda,P);
    V = P - Z;
end

% horizontal processing
W = rot90(X);
for k = 1:N
    X = M\(W+gamma\Psi'*(Z-V));
    P = Psi*X+V;
    Z = soft_thresholding(gamma*lambda,P);
    V = P - Z;
end
X = rot90(X,3);

%% Result
figure;
imshow(uint8(round(X))); % for lena
title('Restored image');
