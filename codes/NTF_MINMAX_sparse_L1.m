function [ntf, Rf, gopt, diagn] = NTF_MINMAX_sparse_L1(lambda, order, Omega, H_inf, zf, threshold)
% [ntf, Rf, gopt, diagn] = NTF_MINMAX_sparse_L1(lambda, order, OSR, H_inf, f0, zf)
%Synthesize a noise transfer function (NTF) for a lowpass or bandpass delta-sigma modulator
%by min-max optimization.
% Arguments:
%   lambda: L1 regularization parameter
%       The default value is lambda=1.
%	order:  The order of NTF which is an FIR filter. 
%		The default value is order=32.
%	omega:  Cut-off frequency
%		The default value is pi/32;
%	H_inf:  The maximum out-of-band gain of NTF. 
%		The default value is Hinf=1.5.
% 	zf:     The flag for assigning NTF zeros. 
%		If zf is 0, then the design is executed without zero assignment. 
%		Otherwise, zeros of the NTF to be optimized is assigned at the center frequency. 
%		The default value is zf=0.
%   threshold:  sparse vector threshold
%
% OUTPUTS:
%	ntf:	The optimized NTF given as a zpk object. See zpk.m
%	Rf:	FIR filter coefficients of the optimized loop filter R(z)
% NOTICE:
%	This function uses free toolboxes YALMIP and SDPT3
%	For YALMIP, see the following web page:
%		https://yalmip.github.io/
%	For SDPT3, see
%		https://github.com/SQLP/SDPT3
%

% State space representation of F
A = [zeros(order-1,1) , eye(order-1) ; zeros(1,order)];
B = [zeros(order-1,1) ; 1];

% Definition of LMI's and LME's
    % LMI for gain optimization
    c = sdpvar(1,order);
    P = sdpvar(order, order, 'symmetric');
    Q = sdpvar(order, order, 'symmetric');
    g = sdpvar(1,1);
    M1 = A'*P*A+Q*A+A'*Q-P-2*Q*cos(Omega);
    M2 = A'*P*B + Q*B;
    M3 = B'*P*B-g;
    M = [M1,M2,c';M2',M3,1;c,1,-1];
    Constraints = [Q>=0, M<=0];
    if zf==1 %Placing a zero at z=1
        Constraints = [Constraints, sum(c)==-1];
    end
    if H_inf<inf %H-infinity norm condition of NTF for stability
        R = sdpvar(order,order,'symmetric');
        Constraints = [Constraints, R>=0, [A'*R*A-R, A'*R*B, c';B'*R*A, -H_inf^2+B'*R*B,1;c,1,-1]<=0];
    end

% Solve LMI's and LME's
Constraints = [Constraints, g>=0];
Objective = g + lambda * norm(c,1);

diagn=optimize(Constraints,Objective,sdpsettings('solver','sdpt3'));

% thresholding
c_vec = double(c);
c_vec(abs(c_vec)<threshold) = 0;

% Outputs
Rf = [0,fliplr(c_vec)]; % Filter coefficients
NTFss = ss(A,B,c_vec,1,1);
ntf = zpk(NTFss); % NTF
gopt = sqrt(double(g)); % Optimal value

return;