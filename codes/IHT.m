
function [x,nitr]=IHT(y,Phi,lambda,gamma,EPS,MAX_ITER)
%IHT iterative hard thresholding algorithm

    if nargin==2
        lambda = 1;
        gamma = 1;
        EPS = 1e-6;
        MAX_ITER = 1000;
    elseif nargin==3
        gamma = 1;
        EPS = 1e-6;
        MAX_ITER = 1000;
    elseif nargin==4
        EPS = 1e-6;
        MAX_ITER = 1000;
    end
    
    [m,n] = size(Phi);
    x = zeros(n,1);
    r = y;
    k = 0;
    while (norm(r)>EPS) & (k < MAX_ITER)
        p = x + gamma * Phi'*r;
        x = hard_thresholding(sqrt(2*lambda*gamma),p);
        S = supp(x);
        r = y-Phi(:,S)*x(S);%norm(r)
        k = k+1;
    end
    nitr=k;
end