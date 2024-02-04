clear
close all

tspan = [0 1500];

a = @(t) t./600;

odefun = @(t, y) [5*y(1)*(1-y(1))*(y(1)-0.1)-y(2)+a(t);...
    0.1*(y(1)-0.25*y(2))];

y0 = [0 0]';

[t, u] = ode15s(odefun, tspan, y0);

figure(1)
subplot(1,2,1)
plot(t, u)

subplot(1,2,2)
plot(u(:,1), u(:,2))