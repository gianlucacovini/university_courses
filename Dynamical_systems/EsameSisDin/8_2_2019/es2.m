clear
close all

tspan = [0 25];

y0 = [3 3 1]';

A = [-1 20 0; -20 -1 0; 0 0 0.2];
odefun = @(t, y) A*y;

%% Simulazione con ode45

[~, y45] = ode45(odefun, tspan, y0);

figure(1)
plot3(y45(:,1), y45(:,2), y45(:,3))
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

%% Simulazione con ode113

[~, y113] = ode113(odefun, tspan, y0);

figure(2)
plot3(y113(:,1), y113(:,2), y113(:,3))
title('Simulazione con ode113')
xlabel('y1')
ylabel('y2')
zlabel('y3')

%% Dato iniziale (-3, -3, -1)

y0 = [-3 -3 -1]';

%% Simulazione con ode45

[~, y45] = ode45(odefun, tspan, y0);

figure(3)
plot3(y45(:,1), y45(:,2), y45(:,3))
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

%% Simulazione con ode113

[~, y113] = ode113(odefun, tspan, y0);

figure(4)
plot3(y113(:,1), y113(:,2), y113(:,3))
title('Simulazione con ode113')
xlabel('y1')
ylabel('y2')
zlabel('y3')