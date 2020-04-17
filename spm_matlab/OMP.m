
function [x,nitr]=OMP(y,Phi,EPS,MAX_ITER)

    if nargin==2
        EPS = 1e-6;
        MAX_ITER = 1000;
    elseif nargin==3
        MAX_ITER = 1000;
    end
    
    Phi_norm = diag(Phi'*Phi);
    
    [m,n] = size(Phi);
    x = zeros(n,1);
    r = y;
    k = 0;
    S = zeros(n,1);
    while  (norm(r)>EPS) & (k < MAX_ITER)
        p = Phi'*r;
        v = p./sqrt(Phi_norm);
        [z,ik] = max(abs(v));
        S(ik) = ik;
        Phi_S = Phi(:,S(S>0));
        %x(S(S>0)) = Phi_S\y;
        x(S(S>0)) = pinv(Phi_S)*y;
        r = y-Phi*x;
        k = k+1;
    end
    nitr=k;
end