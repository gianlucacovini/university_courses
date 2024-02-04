function u = secondoMetodo(h, f, y0, N, alpha)
    u = zeros(N, 2);
    u(1, :) = y0;
    for n=1:N-1
        k1 = f(u(n, 1), u(n, 2));
        k2 = f(u(n,1)+alpha*h(n)*k1(1), u(n,2)+alpha*h(n)*k1(2));
        u(n+1, :) = u(n, :) + h(n)*((1-1/(2*alpha))*k1+1/(2*alpha)*k2);
    end
end