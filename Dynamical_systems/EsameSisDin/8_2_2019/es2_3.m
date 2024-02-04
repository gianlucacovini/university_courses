clear
close all

tspan = [0 25];

y0 = [3 3 1]';

A = [-1 20 0; -20 -1 0; 0 0 -0.2];
odefun = @(t, y) A*y;

%% Simulazione con ode45

[~, y45pos] = ode45(odefun, tspan, y0);

%% Simulazione con ode113

[~, y113pos] = ode113(odefun, tspan, y0);

%% Dato iniziale (-3, -3, -1)

y0 = [-3 -3 -1]';

%% Simulazione con ode45

[~, y45neg] = ode45(odefun, tspan, y0);

%% Simulazione con ode113

[~, y113neg] = ode113(odefun, tspan, y0);

%% Grafici

figure(1)
subplot(2,2,1)
plot3(y45pos(:,1), y45pos(:,2), y45pos(:,3))
grid on
grid minor
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,2)
plot3(y113pos(:,1), y113pos(:,2), y113pos(:,3))
grid on
grid minor
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,3)
plot3(y113neg(:,1), y113neg(:,2), y113neg(:,3))
grid on
grid minor
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')

subplot(2,2,4)
plot3(y45neg(:,1), y45neg(:,2), y45neg(:,3))
title('Simulazione con ode45')
xlabel('y1')
ylabel('y2')
zlabel('y3')