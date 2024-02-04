clear
close all

N = 10000;
T = 20;
t0 = 0;

y10 = 1;
y20 = -3;
y0 = [y10 y20]';

h0 = (T-t0)/N;
h = ones(1,N)'*h0;

mu = 1;

f = @(y1, y2) [mu*y1+y2-(y1+y2)*(y1^2+y2^2);...
    -y1+mu*y2-(y2-y1)*(y1^2+y2^2)]';

u = primoMetodo(h, f, y0, N);

t = linspace(t0, T, N);
%plot3(t, u(:,1), u(:,2))
plot(u(:,1), u(:,2))

hold on

%% Soluzione esatta

r0 = sqrt((y0(1)^2+y0(2)^2));

theta0 = atan(y0(2)/y0(1));

phi = @(t) r0^2+(mu-r0^2)*exp(-2*mu*t);

r = @(t) (sqrt(mu)*r0)./(sqrt(r0^2+(mu-r0^2)*exp(-2*mu*t)));

theta = @(t) theta0-t+mu*(t+(log(phi(t)-log(mu)))/2*mu);

y = @(t) [r(t).*cos(theta(t)); r(t).*sin(theta(t))];

ysol = y(t);

%plot3(t, ysol(1, :), ysol(2, :))
plot(ysol(1, :), ysol(2, :))

hold on

%% Secondo metodo

alpha = 0.1;

v = secondoMetodo(h, f, y0, N, alpha);

%plot3(t, v(:,1), v(:,2))
plot(v(:,1), v(:,2))