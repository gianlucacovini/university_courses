function u = primoMetodo(h, f, y0, N)
    u = zeros(N, 2);
    u(1, :) = y0;
    for n=1:N-1
        K_1n = f(t(n), u(n, :));
        u(n+1, :) = u(n, :) + h(n)*K_1n;
    end
end