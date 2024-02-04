clear
close all

f = @(t, y) -y;
tspan = [0 400];
y0 = 1;

eps = 0.01;

y = @(t) y0*exp(1).^(-t);

figure(1)
h = 0;
errInf = 0;
while errInf < 1e3
    h = h+eps;
    [u, t] = AB3(f, tspan, y0, h);

    plot(t, u)
    title('Soluzione esatta e approssimata')
    xlabel('t')
    ylabel('y(t)')
    hold on
    
    errInf = norm(u-y(t), 'inf');
end

hmin = h
plot(t, y(t), '--', 'LineWidth', 2)
hold off