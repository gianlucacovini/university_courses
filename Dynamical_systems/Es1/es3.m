clear
close all

tspan = [0 300];
sigma = 10;
r = 28;
b = 8/3;

y10 = 0;
y20 = 10;
y30 = 0;

[t,y] = ode15s(@(t,y) [sigma*(y(2)-y(1)); r*y(1)-y(2)-y(1)*y(3);...
    y(1)*y(2)-b*y(3)], tspan, [y10;y20;y30]);

figure
subplot(4,1,1)
plot(t,y(:,1))
title('Grafico y_1')
xlabel('t')
ylabel('y_1(t)')

subplot(4,1,2)
plot(t,y(:,2))
title('Grafico y_1')
xlabel('t')
ylabel('y_1(t)')


subplot(4,1,3)
plot(t,y(:,3))
title('Grafico y_1')
xlabel('t')
ylabel('y_1(t)')

subplot(4,1,4)
plot(t,sqrt(y(:,1).^2+y(:,2).^2+y(:,3).^2))
title('Grafico r(t)')
xlabel('t')
ylabel('r(t)')

figure

plot3(y(:,1),y(:,2),y(:,3))
title('Orbita (y_1, y_2, y_3)')
xlabel('y_1(t)')
ylabel('y_2(t)')
zlabel('y_3(t)')