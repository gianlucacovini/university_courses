function [u, t] = EE(f, tspan, y0, h)
    t0 = tspan(1);
    T = tspan(2);
    t(1) = t0;
    u(1) = y0;
    for n=1:floor((T-t0)/h)
        u(n+1) = u(n)+h*f(t(n),u(n));
        t(n+1) = t(n)+h;
    end
end