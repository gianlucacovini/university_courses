% Errore calore

%% Poisson
close all
clear 

xspan = [0 1];
alpha = 1;
beta = -1;

num_h = 5;
E = zeros(num_h, 1);
h_vec = zeros(num_h, 1);
c = 1;
for exponent=1:num_h
    n=2^(8 - exponent);
    
    h_vec(c) = 1/(n+1);
    h = 1/(n+1);
    
    A = zeros(n, n);
    for i=1:n
        for j=1:n
            if i == j
                A(i, j) = 2/h^2;
            elseif i == j+1 || i == j-1
                A(i, j) = -1/h^2;
            end
        end
    end
    
    f = ones(n, 1);
    f(1) = f(1) + alpha/h^2;
    f(n) = f(n) + beta/h^2;
    
    U = A\f;
    
    U_tilde = [alpha; U; beta];
    
    points = zeros(n+1, 1);
    for i=1:(n+1)
        points(i) = i*h;
    end
    
    points = [0; points];
    
    u = @(x) -1/2 * x.^2 -3/2*x + 1;
    u_e = u(points);

    E(exponent) = norm(u_e-U_tilde, 'inf');
    
    %figure(1)
    %plot(points, U_tilde, '--', points, u_e)

    c = c+1;
end

figure(2)
loglog(h_vec, E, 'o-', h_vec, h_vec.^2, '--', h_vec, h_vec, '--')
legend('Err', 'h^2', 'h')

%% Poisson (b)

xspan = [0 1];
alpha = 1;
beta = -1;

num_h = 4;
E = zeros(num_h, 1);
h_vec = zeros(num_h, 1);
c = 1;
for exponent=1:num_h
    n=2^(8 - exponent);
    
    h_vec(c) = 1/(n+1);
    h = 1/(n+1);
    
    A = zeros(n, n);
    for i=1:n
        for j=1:n
            if i == j
                A(i, j) = 2/h^2;
            elseif i == j+1 || i == j-1
                A(i, j) = -1/h^2;
            end
        end
    end

    
    points = zeros(n+1, 1);
    for i=1:(n+1)
        points(i) = i*h;
    end
    
    points = [0; points];
    
    f = sin(points(2:end-1));
    f(1) = f(1) + alpha/h^2;
    f(n) = f(n) + beta/h^2;
    
    U = A\f;
    
    U_tilde = [alpha; U; beta];
    
    u = @(x) (-sin(1)-2)*x+1+sin(x);
    u_e = u(points);

    E(exponent) = norm(u_e-U_tilde, 'inf');
    
    %figure(3)
    %plot(points, U_tilde, '--', points, u_e)

    c = c+1;
end

figure(4)
loglog(h_vec, E, 'o-', h_vec, h_vec.^2, '--', h_vec, h_vec, '--')
legend('Err', 'h^2', 'h')

%% Calore

xspan = [0 1];
T_f = 0.2;
tspan = [0 T_f];
alpha = 0;
beta = 0;
n = 100;
m = 50;

h = 1/(n+1);

A = zeros(n, n);
for i=1:n
    for j=1:n
        if i == j
            A(i, j) = 2/h^2;
        elseif i == j+1 || i == j-1
            A(i, j) = -1/h^2;
        end
    end
end

points = zeros(n+1, 1);
for i=1:(n+1)
    points(i) = i*h;
end

delta_t = T_f/m;
times = zeros(m, 1);
for i=1:m
    times(i) = tspan(1)+i*delta_t;
end

U = zeros(m, n);
U(1, :) = 100*sin(pi*points(1:end-1)); 
delta_t = T_f/m;
for k=1:m-1
    I = eye(n);

    U(k+1, :) = (I + delta_t*A)\U(k, :)';
end

[X, T] = meshgrid(points(1:end-1), times);  % Crea una griglia di punti per x e t

figure(5)
surf(X, T, U)  % Plotta i dati come una superficie
xlabel('x')
ylabel('t')
zlabel('u')

%% Errore Calore 
clear

xspan = [0 1];
T_f = 0.2;
tspan = [0 T_f];
alpha = 0;
beta = 0;
m = 100;

num_h = 4;
E = zeros(num_h, 1);
h_vec = zeros(num_h, 1);
c = 1;
for exponent=1:num_h
    n=2^(9 - exponent);
    
    h_vec(c) = 1/(n+1);
    h = 1/(n+1);
    
    A = zeros(n, n);
    for i=1:n
        for j=1:n
            if i == j
                A(i, j) = 2/h^2;
            elseif i == j+1 || i == j-1
                A(i, j) = -1/h^2;
            end
        end
    end
    
    points = zeros(1, n);
    for i=1:n
        points(i) = i*h;
    end

    points = [xspan(1) points xspan(2)];
    
    U = zeros(m, n);
    U(1, :) = 100*sin(pi*points(2:end-1));
    delta_t = T_f/m;
    for k=1:m
        I = eye(n);
    
        U(k+1, :) = (I + delta_t*A)\U(k, :)';
    end

    % Aggiungere termini di bordo a U (tipo U_tild di prima)  !!
    U = [zeros(m+1, 1) U zeros(m+1, 1)];
    
    times = zeros(m, 1);
    for i=1:m
        times(i) = tspan(1)+i*delta_t;
    end

    times = [tspan(1); times];
    
    e = exp(1);
    u = @(x, t) 100*e.^(-pi^2*t).*sin(pi*x);
    u_e = u(points, times);
    
    for k=1:m
        E(exponent, k) = norm(u_e(k, :)-U(k, :), 'inf');
    end

    c = c+1;
end

[X, T] = meshgrid(points, times);  % Crea una griglia di punti per x e t

figure(6)
surf(X, T, U)  % Plotta i dati come una superficie
xlabel('x')
ylabel('t')
zlabel('u')

[rows_E, ~] = size(E);
E_inf = zeros(1, rows_E);
for i=1:rows_E
    E_inf(i) = norm(E(i, :));
end

figure(7)
loglog(h_vec, E_inf, h_vec, h_vec, '--', h_vec, h_vec.^2, '--')

%%
close all

figure(1)
surf(X, T, u_e)
xlabel('x')
ylabel('t')
zlabel('u')

figure(2)
surf(X, T, U)
xlabel('x')
ylabel('t')
zlabel('u')

figure(3)
surf(X, T, abs(u_e-U))
xlabel('x')
ylabel('t')
zlabel('E')

%% Calore con Neumann

xspan = [0 1];
T_f = 0.2;
tspan = [0 T_f];
alpha = 0;
beta = 0;
n = 100;
m = 50;

h = (xspan(2)-xspan(1))/(n+1);

A = zeros(n+2, n+2);
for i=1:n+2
    for j=1:n+2
        if i == j
            A(i, j) = 2/h^2;
        elseif i == j+1 || i == j-1
            A(i, j) = -1/h^2;
        end
    end
end

A(1, 2) = -2/h^2;
A(n+2, n+1) = -2/h^2;

points = zeros(n+2, 1);
for i=1:(n+2)
    points(i) = xspan(1) + i*h;
end

U = zeros(m, n+2);
u0 = @(x) (x >= 0.4).*(x <= 0.6);
U(1, :) = u0(points);

delta_t = T_f/m;
for k=1:m-1
    I = eye(n+2);

    U(k+1, :) = (I + delta_t*A)\U(k, :)';
end

[X, T] = meshgrid(points(2:end-1), times);  % Crea una griglia di punti per x e t

figure(5)
surf(T, X, U(:, 2:end-1))  % Plotta i dati come una superficie
xlabel('t')
ylabel('x')
zlabel('u')