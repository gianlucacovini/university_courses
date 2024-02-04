function u = quartoMetodo(h, f, y0, N)
    u = zeros(N, 2);
    u(1, :) = y0;
    for n=1:N-1
        k1 = f(u(n, 1), u(n, 2));
        k2 = f(u(n,1)+h(n)*k1(1)/2, u(n,2)+h(n)*k1(2)/2);
        k3 = f(u(n,1)+h(n)*k2(1)/2, u(n,2)+h(n)*k2(2)/2);
        k4 = f(u(n,1)+h(n)*k3(1), u(n,2)+h(n)*k3(2));        
        u(n+1, :) = u(n, :) + h(n)/3*(k1/2+k2+k3+k4/2);
    end
end