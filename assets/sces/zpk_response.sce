// Section 3.4
// Poles and impulse response
// For Scilab
// exec("sec3_4.sci",-1)

clear;

// Define 's'
s = poly(0,'s');
// Define j=sqrt(-1)
j = %i;

//// Poles
// Complex poles
p1 = -1 + j; // j = sqrt(-1)
p2 = -2-2*j;
p1c = conj(p1); // complex conjugate of p1
p2c = conj(p2); // complex conjugate of p2

// Real poles
p3 = -3;

// The poles
p = [p1; p1c; p2; p2c; p3];

disp('The poles are')
disp(p)

px=real(p);
py=imag(p);

//// Zeros
// Complex zeros
z1 = [];
z2 = []; // [] means no more complex zeros
z1c = conj(z1);
z2c = conj(z2);

// Real zeros
z3 = -1;

// The zeros
z = [z1; z1c; z2; z2c; z3];

disp('The zeros are')
disp(z)

zx=real(z);
zy=imag(z);

// Plot the poles and zeros on the complex plane
scf();
rp = norm(p,'inf');
rz = norm(z,'inf');
r = max([rp,rz]);
plot2d(px,py,style=-2,rect=[-2*r,-2*r,2*r,2*r]);
plot2d(zx,zy,style=-9,rect=[-2*r,-2*r,2*r,2*r]);
xtitle('x: poles,  o: zeros','real','imag');
xgrid();

// System with the poles p1 and p2
num=poly(z,'s');
num_c=coeff(num);
den=poly(p,'s');
den_c=coeff(den);
K=num_c(1)/den_c(1);
G=num/den/K;
Gs=syslin("c",G);

// Time responses
Tf=10;
dt=Tf/1000;
t = 0:dt:Tf; //simulation time

// Impulse response
sysi=csim('impulse', t, Gs);
scf(); plot2d(t,sysi,style=2);

// Step response
syss=csim('step', t, Gs);
plot2d(t,syss,style=5);

// Figure Properties
A=gca();
P=A.children.children;
P(1).thickness=3;
P(2).thickness=3;
xgrid();
legend('impulse response','step response',1);

return;
