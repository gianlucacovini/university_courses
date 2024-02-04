clear
close all

% y1 = 0.1e^mu*t
% y2 = e^(-t)[9/10 + 1/10 e^(1+mu)t

mu = -5;
T = 10;
t0 = 0;
N = 1000;

A = [mu 0; 1+mu -1];

y0 = [0.1 1];

h0 = (T-t0)/N;
h = ones(N,1)*h0;

f = @(y) A*y;

u = quartoMetodoVett(h, f, y0, N);

e = exp(1);
y =@(t) [0.1*e.^(mu*t); e.^(-t).*(9/10 + 1/10*e.^((1+mu)*t))];

t = linspace(t0, T, N);
ysol = y(t)';
figure(1)
plot(t, u(:,1), 'b', 'LineWidth', 2)
hold on
plot(t, ysol(:,1), 'y--', 'LineWidth', 2)
hold off
legend('Soluzione approssimata','Soluzione esatta')
xlabel('t')
ylabel('y(t)')
title('Grafico della funzione a passo costante')

%% Punto 2

N = 1000;

for n=1:floor(4/abs(mu)/0.1)
    h(n) = 0.1;
end

for n=floor(4/abs(mu)/0.1)+1:N
    h(n) = (T-4/abs(mu))/(N-floor(4/abs(mu)/0.1));
end

v = quartoMetodoVett(h, f, y0, N);

t = linspace(t0, T, N);
ysol = y(t)';
figure(2)
plot(t, u(:,1), 'r', 'LineWidth', 2)
hold on
plot(t, ysol(:,1), 'g--', 'LineWidth', 2)
hold off
legend('Soluzione approssimata','Soluzione esatta')
xlabel('t')
ylabel('y(t)')
title('Grafico della funzione a passo variabile')

%% Pto 3

figure(3)
for N=1:50
    for n=1:floor(4/abs(mu)/0.1)
        h(n) = 0.1;
    end

    for n=floor(4/abs(mu)/0.1):N
        h(n) = (T-4/abs(mu))/(N-floor(4/abs(mu)/0.1));
    end

    v = quartoMetodoVett(h, f, y0, N);

    t = linspace(t0, T, N);
    ysol = y(t)';
    semilogy(t, v(:,1))
    hold on
end

hold off

% N = 25, h = 10-0.8/25 = 0.368

%% pto 4

N = 100;

for n=1:floor(4/abs(mu)/0.1)
    h(n) = 0.1;
end

for n=floor(4/abs(mu)/0.1):N
    h(n) = (T-4/abs(mu))/(N-floor(4/abs(mu)/0.1));
end

v = quartoMetodoVett(h, f, y0, N);

t = linspace(t0, T, N);
ysol = y(t)';

err1 = abs(v(:,1)-ysol(:,1));

err2 = abs(v(:,2)-ysol(:,2));

Einf = zeros(1,N);
E2 = zeros(1,N);

for n=1:N
    Einf(n) = max(abs([err1(n),err2(n)]));
    E2(n) = sqrt(err1(n)^2+err2(n)^2);
end

%% Grafici

figure(4)
subplot(2,2,1)
plot(v(:,1), ysol(:,1),'o', 'LineWidth', 2)
hold on
plot(v(:,2), ysol(:,2), 'LineWidth', 2)
hold off

subplot(2,2,2)
plot(t,err1,'r-o', 'LineWidth', 2)
hold on
plot(t,err2,'b-o', 'LineWidth', 2)
hold off

norm2e = zeros(1,N);
norm2a = zeros(1,N);
for n=1:N
    norm2e(n) = norm(ysol(n,:));
    norm2a(n) = norm(v(n,:));
end

subplot(2,2,3)
plot(t,norm2e, 'g-o', 'LineWidth', 2)
hold on
plot(t,norm2a, 'LineWidth', 2)

subplot(2,2,4)
plot(ysol(:,1),ysol(:,2), 'g', 'LineWidth', 2)
hold on
plot(v(:,1),v(:,2), 'b--', 'LineWidth', 2)
hold off     