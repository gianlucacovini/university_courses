clear
close all

f = @(t, y) -100*y+100*t+1;
tspan = [0 10];
y0 = 3;

y = @(t) y0*exp(1).^(-100*t)+t;

%% EE h=0.015

h = 0.015;

[uE1, tE1] = EE(f, tspan, y0, h);

numItE1 = size(tE1, 2);
for n=1:numItE1
    EE1(n) = uE1(n)-y(tE1(n));
end

%% EE h=0.03

h = 0.03;

[uE2, tE2] = EE(f, tspan, y0, h);

numItE2 = size(tE2, 2);
for n=1:numItE2
    EE2(n) = uE2(n)-y(tE2(n));
end

%% Grafici

figure(1)
subplot(2,1,1)
plot(tE1, uE1, tE1, y(tE1), '--', 'LineWidth', 2)
title('Eulero esplicito con h=0.015')
legend('Soluzione approssimata', 'Soluzione esatta', 'Location', 'northwest')

subplot(2,1,2)
plot(tE2, uE2, tE2, y(tE2), '--', 'LineWidth', 2)
title('Eulero esplicito con h=0.03')
legend('Soluzione approssimata', 'Soluzione esatta', 'Location', 'northwest')