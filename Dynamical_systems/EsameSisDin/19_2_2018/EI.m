function [u, t] = EI(f, df, tspan, y0, h)
    t0 = tspan(1);
    T = tspan(2);
    nsteps = floor((T-t0)/h);
    maxiter = 20;
    tol = 1e-10;
    B = [1 1; 0 1];
    C = B(1:end-1,1);
    b = B(end, 2:end)';
    A = B(1:end-1,2:end);
    t = linspace(t0, T, nsteps);
    Y = zeros(nsteps,1);
    h = t(2)-t(1);
    Y(1) = y0;
    nstages = size(B,1)-1;
    for n = 1:nsteps-1 
        Kold = zeros(nstages,1);
        iter = 0;
        err = 1;
        while (err> tol) && (iter<maxiter)
            iter = iter+1;
            num = Kold - f(t(n),Y(n)+h*A*Kold);
            DF = 1-h*A*df(t(n),Y(n)+h*A*Kold);
            dif = DF\num;
            Knew = Kold - dif;
            Kold = Knew;
            err = max(abs(dif));
        end
        it(n) = iter;
        Y(n+1) = Y(n) + h*(b'*Knew);
    end
    u = Y;
end