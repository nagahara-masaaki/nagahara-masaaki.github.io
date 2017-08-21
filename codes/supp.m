
function I = supp(x)
%SUPP support of a vector
% SUPP(X) returns the support set (the index set on which X is non-zero).
    I = find(abs(x)>0)';
end