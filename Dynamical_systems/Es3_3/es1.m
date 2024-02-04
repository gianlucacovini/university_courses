clear
close all

for i=1:4
    f = @(t,y) ones(size(y)).*y-t+1;
    df = @(t,y) ones(size(y)).*1;
    y0 = 1;
    T = 1;
    nsteps = 100*2^i;
    maxiter = 20;
    tol = 1e-10;
    B = [1 1; 0 1];
    C = B(1:end-1,1);
    b = B(end, 2:end)';
    A = B(1:end-1,2:end);
    t = linspace(0, T, nsteps);
    Y = zeros(nsteps,1);
    h = t(2)-t(1);
    H(i) = i;
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
    e = exp(1);
    yex = @(t) e.^t+t;
   
    E(i) = norm(yex(t)-Y')/norm(yex(t));
    I(i) = sum(it)/n;
end

subplot(2,2,1)
plot(t, Y, 'o-',t, yex(t), 'r')

subplot(2,2,2)
plot(H,E)

subplot(2,2,3)
plot(H,I)