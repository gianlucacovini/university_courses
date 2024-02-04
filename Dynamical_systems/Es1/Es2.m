clear
close all

tspan = 0:0.01:100;
u = @(t) exp(1).^(-t)+sin(t)+cos(sqrt(2).*t);

figure
plot(tspan,u(tspan))
title('Funzione u')
xlabel('t')
ylabel('u(t)')