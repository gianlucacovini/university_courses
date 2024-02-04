function u = primoMetodo(h, f, y0, N)
    u = zeros(N, 2);
    u(1, :) = y0;
    for n=1:N-1
        u(n+1, :) = u(n, :) + h(n)*f(u(n, 1), u(n, 2));
    end
end