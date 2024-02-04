clear
close all

tspan = [0 20];

y0 = [3 3 1]';

A = [-1 20 0; -20 -1 0; 0 0 -0.2];
odefun = @(t, y) A*y;

%% Simulazione con ode45

[t45, y45] = ode45(odefun, tspan, y0);
NumIt45 = size(t45, 1);

h45 = zeros(NumIt45-1, 1);
for i=1:NumIt45-1
    h45(i) = t45(i+1)-t45(i);
end

hmin45 = min(h45);
hmax45 = max(h45);

T45 = t45(1:NumIt45-1);

%% Simulazione con ode113

[t113, y113] = ode113(odefun, tspan, y0);
NumIt113 = size(t113, 1);

h113 = zeros(NumIt113-1, 1);
for i=1:NumIt113-1
    h113(i) = t113(i+1)-t113(i);
end

hmin113 = min(h113);
hmax113 = max(h113);

T113 = t113(1:NumIt113-1);

%% Grafici

figure(1)
subplot(2,1,1)
plot(T45, h45)
title('ode45')
xlabel('t')
ylabel('h')

subplot(2,1,2)
plot(T113, h113)
title('ode113')
xlabel('t')
ylabel('h')