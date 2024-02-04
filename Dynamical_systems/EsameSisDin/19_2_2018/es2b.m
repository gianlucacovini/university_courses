clear
close all

tspan = [0 100];

b = 1;
y0 = 0.99;
f = @(t, y) -y*sqrt(b-y);

df = @(t, y) (3*y-2)/(2*sqrt(1-y));

h = 0.01;
[tE, uE] = EE(f, tspan, y0, h);

[tI, uI] = EI(f, df, tspan, y0, h);

figure(1)
subplot(2,1,1)
plot(tE, uE)

subplot(2,1,2)
plot(tI, uI)
