clear
close all

tspan = [0 1];
y0 = 0;

h = 1/100;

e = exp(1);
f = @(t, y) -t*e^(-y);

u(1) = y0;
t = 1;
for n=1:100
    u(n+1) = u(n)+h*f(t,u(n));
end

u100 = u(100);