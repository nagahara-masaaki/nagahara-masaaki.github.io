
function [x,nitr]=MP(y,Phi,EPS,MAX_ITER)

    if nargin==2
        EPS = 1e-6;
        MAX_ITER = 1000;
    elseif nargin==3
        MAX_ITER = 1000;
    end
    
    
    [m,n] = size(Phi);
    x = zeros(n,1);
    r = y;
    k = 0;
    Phi_norm = diag(Phi'*Phi);
    while (norm(r)>EPS) & (k < MAX_ITER)
        p = Phi'*r;
        v = p./sqrt(Phi_norm);
        [z,ik] = max(abs(v));
        v2 = p./Phi_norm;
        z = v2(ik);
        x(ik) = x(ik)+z;
        r = r-z*Phi(:,ik);
        k = k+1;
    end
    nitr=k;
end