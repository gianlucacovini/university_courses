clear
close all

% ro(x) = x - 1;
% sigma(x) = 1/2 x +1/2;

f = @(Z) (Z-1)./(1/2 *Z+1/2);

theta = linspace(0, 2*pi, 1000);
Z = exp(1).^(1i*theta);
F = f(Z);

figure(1)
plot(F)
grid on
grid minor
title('Regione di assoluta stabilità')
xlabel('Re')
ylabel('Im')