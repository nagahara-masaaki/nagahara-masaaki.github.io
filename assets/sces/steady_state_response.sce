// Steady-state response
clear;

// Define 's'
s = poly(0,'s');

// Plant P(s)
P1 = 1/(s+1); //for Example 1
P2 = 1/s; // for Example 2

// Controller K(s)
// for Example 1
K1 = 1;
K2 = 3;
K3 = 9;

// for Example 2
K4 = 5/(s+1);
K5 = 1/(s+1);
K6 = 0.5/(s+1);

// Feedback system
// for Example 1
G1 = syslin("c", P1*K1/(1+P1*K1));
G2 = syslin("c", P1*K2/(1+P1*K2));
G3 = syslin("c", P1*K3/(1+P1*K3));

// for Example 2
G4 = syslin("c", P2*K4/(1+P2*K4));
G5 = syslin("c", P2*K5/(1+P2*K5));
G6 = syslin("c", P2*K6/(1+P2*K6));

// Step response
// Feedback system
// Example 1
t1 = 0:0.01:5;
G3step = csim('step', t1, G3);
G2step = csim('step', t1, G2);
G1step = csim('step', t1, G1);
scf();
plot2d(t1,[G3step', G2step', G1step']);
xtitle('Example 4.2; P(s) = 1/(s+1), K(s) = K','t','y(t)');
legend('K=9','K=3','K=1',4);

//Example 2
t2 = 0:0.01:15;
G4step = csim('step', t2, G4);
G5step = csim('step', t2, G5);
G6step = csim('step', t2, G6);
scf();
plot2d(t2,[G4step',G5step',G6step']);
legend('K_0=5','K_0=1','K_0=0.5',4);
xtitle('Example 4.3; P(s) = 1/s, K(s) = K_0/(s+1)','t','y(t)');