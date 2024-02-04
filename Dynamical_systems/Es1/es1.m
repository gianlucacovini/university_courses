clear
close all

mu = 1;
% mu = 0.5

tspan = [0 20];
y1_0 = 3;
y2_0 = -3;

[t,y] = ode15s(@(t,y) [mu*y(1)-y(2)-y(1)*(y(1)^2+3/2 *y(2)^2);...
    y(1)+mu*y(2)-y(2)*(y(1)^2+1/2*y(2)^2)], tspan, [y1_0;y2_0]);

figure(2)
subplot(2,2,1)
plot(t,y(:,1))
title('Grafico y1')
xlabel('t')
ylabel('y1(t)')

subplot(2,2,2)
plot(t,y(:,2))
title('Grafico y2')
xlabel('t')
ylabel('y2(t)')

subplot(2,2,3)
plot(y(:,2),y(:,1))
title('Orbita y2-y1')
xlabel('y1(t)')
ylabel('y2(t)')

subplot(2,2,4)
plot(t,sqrt(y(:,1).^2+y(:,2).^2))
title('Grafico r(t)')
xlabel('t')
ylabel('r(t)')