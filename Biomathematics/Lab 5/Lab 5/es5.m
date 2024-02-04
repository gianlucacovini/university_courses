% Capire cosa chiede punto a

%% Kolmogorov-Fisher

xspan = [0 1];
T_f = 10;
tspan = [0 T_f];
n = 100;
m = 50;

sigma = 1e-3;
b = 5;
alpha = 0.1;
I_app_fun = @(t, x) (x' >= 0 & x' <= 0.04) & (t >= 0 & t <= 1);
% Gives a matrix with the 1-0 times on the rows and points on the columns

h = (xspan(2)-xspan(1))/(n+1);

A = zeros(n+2, n+2);
for i=1:n+2
    for j=1:n+2
        if i == j
            A(i, j) = 2;
        elseif i == j+1 || i == j-1
            A(i, j) = -1;
        end
    end
end

A(1, 2) = -2;
A(n+2, n+1) = -2;
A = sigma*A/h^2;

points = zeros(n+2, 1);
for i=1:(n+2)
    points(i) = xspan(1) + (i-1)*h;
end

U = zeros(m, n+2);
u0 = @(x) 0;
U(1, :) = u0(points);

delta_t = T_f/m;

times = zeros(m+1, 1);
for i=1:m+1
    times(i) = tspan(1) + delta_t*(i-1);
end

I_app_vec = I_app_fun(times, points);
for k=1:m
    I = eye(n+2);

    U(k+1, :) = (I + delta_t*A)\(U(k, :)'+b*U(k, :)'.*(1-U(k, :)')+I_app_vec(k, :)');
end

[X, T] = meshgrid(points, times);  % Crea una griglia di punti per x e t

figure(1)
surf(T, X, U)  % Plotta i dati come una superficie
xlabel('t')
ylabel('x')
zlabel('u')



