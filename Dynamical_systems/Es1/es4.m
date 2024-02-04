clear
close all



T = 300;
tspan = [0 T];
ro = 0.25;
alpha = pi/4;

y10 = 1-ro;
y20 = 0;
y30 = 0;
y40 = alpha*((1+ro)/(1-ro))^(1/2);

opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);

[t,y] = ode15s(@(t,y) [y(3); y(4); -y(1)/((y(1)^2+y(2)^2)^(3/2)/alpha^2);...
    -y(2)/((y(1)^2+y(2)^2)^(3/2)/alpha^2)], tspan, [y10;y20;y30;y40], opts_1);

figure
subplot(2,2,1)
plot(t,y(:,1))
title('Grafico y_1')
xlabel('t')
ylabel('y_1(t)')

subplot(2,2,2)
plot(t,y(:,2))
title('Grafico y_2')
xlabel('t')
ylabel('y_2(t)')

subplot(2,2,3)
plot(t,sqrt(y(:,1).^2+y(:,2).^2))
title('Grafico r(t)')
xlabel('t')
ylabel('r(t)')

subplot(2,2,4)
plot(y(:,1),y(:,2))
title('Orbita (y_1,y_2)')
xlabel('y_1(t)')
ylabel('y_2(t)')