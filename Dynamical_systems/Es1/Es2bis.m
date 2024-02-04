clear
close all

a = pi/3.1;
% a = sqrt(2);
tspan = 0:0.01:100;
u1 = @(t) exp(1).^(-t)+sin(t);
u2 = @(t) exp(1).^(-t) + cos(a.*t);

figure
plot(u1(tspan),u2(tspan))
title('Orbita u_1 contro u_2')
xlabel('u_1(t)')
ylabel('u_2(t)')