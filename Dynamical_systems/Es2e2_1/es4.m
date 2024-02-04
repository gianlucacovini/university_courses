clear
close all

%% Metodo 1

N = 10000;
T = 20;
t0 = 0;

y10 = 1;
y20 = -3;
y0 = [y10 y20]';

h0 = (T-t0)/N;
h = ones(1,N)'*h0;

t = linspace(t0, T, N);

mu = 1;

f = @(y1, y2) [mu*y1+y2-(y1+y2)*(y1^2+y2^2);...
    -y1+mu*y2-(y2-y1)*(y1^2+y2^2)]';

u = primoMetodo(h, f, y0, N);

%% Metodo 2

alpha = 0.1;

v = secondoMetodo(h, f, y0, N, alpha);

%% Sol esatta

r0 = sqrt((y0(1)^2+y0(2)^2));

theta0 = atan(y0(2)/y0(1));

phi = @(t) r0^2+(mu-r0^2)*exp(-2*mu*t);

r = @(t) (sqrt(mu)*r0)./(sqrt(r0^2+(mu-r0^2)*exp(-2*mu*t)));

theta = @(t) theta0-t+mu*(t+(log(phi(t)-log(mu)))/2*mu);

y = @(t) [r(t).*cos(theta(t)); r(t).*sin(theta(t))];

ysol = y(t);

%% Errori

err1 = u(:,1)-ysol(1,:)';

err2 = u(:,2)-ysol(2,:)';

Einf = zeros(1,N);
E2 = zeros(1,N);

for n=1:N
    Einf(n) = max(abs([err1(n),err2(n)]));
    E2(n) = sqrt(err1(n)^2+err2(n)^2);
end

%% Grafici

subplot(2,2,1)
plot(t, u(:,1),t, ysol(1,:),t, u(:,2),t, ysol(2,:))

subplot(2,2,2)
plot(t,err1,t,err2,t,Einf,t,E2)

norm2e = zeros(1,N);
norm2a = zeros(1,N);
for n=1:N
    norm2e(n) = norm(ysol(:,n));
    norm2a(n) = norm(u(n,:));
end

subplot(2,2,3)
plot(t,norm2e,t,norm2a)

subplot(2,2,4)
plot(u(:,1),u(:,2),ysol(1,:),ysol(2,:))