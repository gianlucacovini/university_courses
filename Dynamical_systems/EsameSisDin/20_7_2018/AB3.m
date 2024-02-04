function [u, t] = AB3(f, tspan, y0, h)
    t0 = tspan(1);
    T = tspan(2);
    t(1) = t0;
    u(1) = y0;

    for n=1:2
        u(n+1) = u(n)+h*f(t(n),u(n));
        t(n+1) = t(n)+h;
    end

    for n=3:floor((T-t0)/h)
        u(n+1) = u(n)+h*(23/12*f(t(n),u(n))-16/12*f(t(n-1),u(n-1))+5/12*f(t(n-2),u(n-2)));
        t(n+1) = t(n)+h;
    end
end