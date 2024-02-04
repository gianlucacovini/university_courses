clear
close all

f = @(t, y) sin(t)+y;
tspan = [0 1];
y0 = 0;

y = @(t) (exp(1).^t-sin(t)-cos(t))./2;

for i=1:10
    nsteps = 2^i;
    [u, t] = metodo(f, tspan, y0, nsteps);

    E(i) = norm(u-y(t))/norm(y(t));
    N(i) = nsteps;
end

dom = 2.^-(1:10);
figure(1)
loglog(dom, E, 'bo-', dom, dom.^3)
title('Grafico convergenza')
xlabel('n')
ylabel('Errore')
legend('Errore in norma 2', 'h^3')