clear
close all

tspan = [0 10];
y0 = [2 0 0];

%% eps=1

eps = 1;

odefun = @(t, y) [(-5*y(1)-y(1)*y(2)+5*y(2)^2+y(3))/eps+y(2)*y(3)-y(1);...
    (10*y(1)-y(1)*y(2)-10*y(2)^2+y(3))/eps-y(2)*y(3)+y(1);...
    (y(1)*y(2)-y(3))/eps-y(2)*y(3)+y(1)];

[t113, y113] = ode113(odefun, tspan, y0);
numIt113 = size(t113, 1);

[t15s, y15s] = ode15s(odefun, tspan, y0);
numIt15s = size(t15s, 1);

%% eps=10^-2

eps = 1e-2;

[t113bis, y113bis] = ode113(odefun, tspan, y0);
numIt113bis = size(t113bis, 1);

[t15sBis, y15sBis] = ode15s(odefun, tspan, y0);
numIt15sBis = size(t15sBis, 1);

% è stiff??

% Come si discute la precisione??

%% Grafici a caso

subplot(2,2,1)
plot(t113, y113)

subplot(2,2,2)
plot(t15s, y15s)

subplot(2,2,3)
plot(t113bis, y113bis)

subplot(2,2,4)
plot(t15sBis, y15sBis)