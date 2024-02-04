clear
close all

%% Adam-Bashford ordine 1

% ro(x) = x - 1;
% sigma(x) = 1;

f = @(Z) Z-1;

theta = linspace(0, 2*pi, 1000);
Z = exp(1).^(1i*theta);
F = f(Z);

figure(1)
plot(F)
grid on
grid minor
title('Adam-Bashfort di ordine 1')
xlabel('Re')
ylabel('Im')

%% Adam-Bashford ordine 2

% ro(x) = x^2 - x;
% sigma(x) = 3/2 x - 1/2;

f = @(Z) (Z.^2-Z)./(3/2*Z-1/2);

theta = linspace(0, 2*pi, 1000);
Z = exp(1).^(1i*theta);
F = f(Z);

figure(2)
plot(F)
grid on
grid minor
title('Adam-Bashfort di ordine 2')
xlabel('Re')
ylabel('Im')

%% Adam-Bashford ordine 3

% ro(x) = x^3 - x^2;
% sigma(x) = 23/12 x^2 - 16/12 x + 5/12;

f = @(Z) (Z.^3-Z.^2)./(23/12*Z.^2-16/12*Z+5/12);

theta = linspace(0, 2*pi, 1000);
Z = exp(1).^(1i*theta);
F = f(Z);

figure(3)
plot(F)
grid on
grid minor
title('Adam-Bashfort di ordine 3')
xlabel('Re')
ylabel('Im')