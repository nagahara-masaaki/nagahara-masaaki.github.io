
function [x,nitr]=iterative_s_sparse(y,Phi,s,gamma,EPS,MAX_ITER)
%ITERATIVE_S_SPARSE iterative s-sparse algorithm

    if nargin==3
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
        x = s_sparse_operator(p,s);
        S = supp(x);
        r = y-Phi(:,S)*x(S);
        k = k+1;
    end
    nitr=k;
end