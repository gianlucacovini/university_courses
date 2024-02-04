clear
close all

T = 20;
t0 = 0;

y10 = 1;
y20 = -3;
y0 = [y10 y20]';

err = zeros(1,5);

for n=1:10
    N = 2^(n+5);
    
    h0 = (T-t0)/N;
    h = ones(1,N)'*h0;

    t = linspace(t0, T, N);

    mu = 1;

    f = @(y1, y2) [mu*y1+y2-(y1+y2)*(y1^2+y2^2);...
        -y1+mu*y2-(y2-y1)*(y1^2+y2^2)]';

    u = primoMetodo(h, f, y0, N);

    alpha = 0.1;

    v = secondoMetodo(h, f, y0, N, alpha);
    
    w = quartoMetodo(h, f, y0, N);

    r0 = sqrt((y0(1)^2+y0(2)^2));

    theta0 = atan(y0(2)/y0(1));

    phi = @(t) r0^2+(mu-r0^2)*exp(-2*mu*t);

    r = @(t) (sqrt(mu)*r0)./(sqrt(r0^2+(mu-r0^2)*exp(-2*mu*t)));

    theta = @(t) theta0-t+mu*(t+(log(phi(t)-log(mu)))/2*mu);

    y = @(t) [r(t).*cos(theta(t)); r(t).*sin(theta(t))];

    ysol = y(t);

    err1 = u(:,1)-ysol(1,:)';
    errb1 = v(:,1)-ysol(1,:)';
    err41 = w(:,1)-ysol(1,:)';    

    err2 = u(:,2)-ysol(2,:)';
    errb2 = v(:,2)-ysol(2,:)';
    err42 = w(:,2)-ysol(2,:)';    

    Einf = zeros(1,N);
    E2 = zeros(1,N);
    Ebinf = zeros(1,N);
    Eb2 = zeros(1,N);

    for i=1:N
        Einf(i) = max(abs([err1(i),err2(i)]));
        E2(i) = sqrt(err1(i)^2+err2(i)^2);
        Ebinf(i) = max(abs([errb1(i),errb2(i)]));
        Eb2(i) = sqrt(errb1(i)^2+errb2(i)^2);   
        E42(i) = sqrt(err41(i)^2+err42(i)^2);        
    end

    norm2e = zeros(1,N);
    norm2a = zeros(1,N);
    normb2a = zeros(1,N);
    for i=1:N
        norm2e(i) = norm(ysol(:,i));
        norm2a(i) = norm(u(i,:));
        normb2a(i) = norm(v(i,:));       
    end
    
    err(n) = norm(E2)/norm(norm2e);
    errb(n) = norm(Eb2)/norm(norm2e);
    err4(n) = norm(E42)/norm(norm2e);    
end

dom = 2.^-(5:14);
loglog(dom, err, dom, dom.*50, dom, dom.^2.*1e4,dom, errb, dom, err4.*10,'bo-', dom, dom.^4.*1e10)
legend('Eulero','h','h^2','R-K due passi','R-K quattro passi','h^4','Location','best')