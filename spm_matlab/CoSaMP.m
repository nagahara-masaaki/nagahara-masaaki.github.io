
function [x,nitr]=CoSaMP(y,Phi,s,EPS,MAX_ITER)

    if nargin==3
        EPS = 1e-6;
        MAX_ITER = 1000;
    elseif nargin==4
        MAX_ITER = 1000;
    end
    
    [m,n] = size(Phi);
    x = zeros(n,1);
    r = y;
    k = 0;
    S = [];
    Lambda = [];
    Phi_norm = diag(Phi'*Phi);
    while (norm(r)>EPS) & (k < MAX_ITER)
        p = s_sparse_operator((Phi'*r)./sqrt(Phi_norm),2*s);
        Ik = supp(p);
        S = union(Lambda,Ik);
        Phi_S = Phi(:,S);
        z = zeros(n,1);
        z(S) = pinv(Phi_S)*y;
        x = s_sparse_operator(z,s);
        Lambda = supp(x);
        r = y-Phi_S*z(S);
        k = k+1;
    end
    nitr=k;
end