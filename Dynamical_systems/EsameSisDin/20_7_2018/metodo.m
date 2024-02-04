function [u, t] = metodo(f, tspan, y0, nsteps)
    t0 = tspan(1);
    T = tspan(2);
    t(1) = t0;
    u(1) = y0;
    h = (T-t0)/nsteps;
    for n=1:nsteps
        k1 = f(t(n), u(n));
        k2 = f(t(n)+1/2*h, u(n)+h/2*k1);
        k3 = f(t(n)+h, u(n)-h*k1+2*h*k2);
        u(n+1) = u(n)+h*(1/6*k1+2/3*k2+1/6*k3);
        t(n+1) = t(n)+h;
    end
end
