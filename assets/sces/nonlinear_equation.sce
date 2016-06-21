//Plot
t=0:0.1:2;
plot(t,t^3-3*t+1)
set(gca(),"grid",[1 1])

//Initial Value
//x0 = 0.5;
x0 = 2;

//Number of Iteration
N = 40;

//Iteration
x = x0;
xa = zeros(1:N+1);
for n = 0 : N
    printf('x[%d] = %f\n',n,x)
    //x = (x^3+1)/3;
    x = (3*x-1)/x^2;
    xa(n+1)=x;
end

figure;
plot(1:N+1,xa);