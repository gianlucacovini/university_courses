clear
close all

tspan = [0 100];

y0 = [1 1 1]';

odefun = @(t, y) [y(2); -y(1); -y(3)];

y = @(t) [sin(t)+cos(t) cos(t)-sin(t) exp(1).^(-t)];

for k=2:9
    options = odeset('RelTol', 10^(-k), 'AbsTol', 10^(-k));

    [t45, u45] = ode45(odefun, tspan, y0, options);
    [t113, u113] = ode113(odefun, tspan, y0, options);
    [t15s, u15s] = ode15s(odefun, tspan, y0, options);
    
    E45(k) = norm(u45-y(t45));
    E113(k) = norm(u113-y(t113));
    E15s(k) = norm(u15s-y(t15s));

    if k==2
        figure(1)
        subplot(2,2,1)
        plot(t45, y(t45), '--')
        hold on
        plot(t45, u45)
        hold off
        title('ode45')

        subplot(2,2,2)
        plot(t113, y(t113), '--')
        hold on
        plot(t45, u45)
        hold off
        title('ode113')

        subplot(2,2,3)
        plot(t15s, y(t15s), '--')
        hold on
        plot(t15s, u15s)
        hold off
        title('ode15s')
    end

    S45(k) = size(t45, 1);
    S113(k) = size(t113, 1);
    S15s(k) = size(t15s, 1);
end

figure(2)
loglog(S45, E45, S113, E113, S15s, E15s, S45, 10^(12)*S45.^(-5), S113, 10^(30)*S113.^(-13))
legend('ode45', 'ode113', 'ode15s', 'h^5', 'h^(13)')
title('Errori in funzioni del numero di passi (es4c)')









