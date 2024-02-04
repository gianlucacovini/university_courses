clear
close all

tspan = [0 25];

y0 = [0 1 0 0.8 0 1.2]';

lambda = @(y) (y(4)^2+y(5)^2+y(6)^2-9.8*y(3))/(y(1)^2+y(2)^2+y(3)^2);
odefun = @(t, y) [y(4); y(5); y(6);...
    -lambda(y)*y(1); -lambda(y)*y(2); -9.8-lambda(y)*y(3)];

options = odeset('RelTol',1e-5);

%% Simulazione con ode23

[t23, y23] = ode23(odefun, tspan, y0, options);
numIt113 = size(t23, 1);

normInf23 = norm(y23(:,1).^2+y23(:,2).^2+y23(:,3).^2+1, 'inf');

%% Simulazione con ode45

[t45, y45] = ode45(odefun, tspan, y0, options);
numIt45 = size(t45, 1);

normInf45 = norm(y45(:,1).^2+y45(:,2).^2+y45(:,3).^2+1, 'inf');

%% Simulazione con ode113

[t113, y113] = ode113(odefun, tspan, y0, options);
numIt113 = size(t113, 1);

normInf113 = norm(y113(:,1).^2+y113(:,2).^2+y113(:,3).^2+1, 'inf');

%% Simulazione con ode15s

[t15s, y15s] = ode15s(odefun, tspan, y0, options);
numIt15s = size(t15s, 1);

normInf15s = norm(y15s(:,1).^2+y15s(:,2).^2+y15s(:,3).^2+1, 'inf');

%% Grafici

figure(1)
subplot(2,2,1)
plot3(y23(:,1), y23(:,2), y23(:,3))
title('Simulazione con ode23')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,2)
plot3(y45(:,1), y45(:,2), y45(:,3))
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,3)
plot3(y113(:,1), y113(:,2), y113(:,3))
title('Simulazione con ode113')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,4)
plot3(y15s(:,1), y15s(:,2), y15s(:,3))
title('Simulazione con ode15s')
xlabel('y1')
ylabel('y2')
zlabel('y3')