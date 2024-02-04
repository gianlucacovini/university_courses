clear
close all

n = 100;
x = linspace(-2.5,2.5,n);
y = x;

%% Primo metodo

[X, Y] = meshgrid(x);

Z = abs(1+X+Y.*1i);

v = [1 1];

figure(1)
contour(X, Y, Z, v)
grid on
grid minor
xlabel('x')
ylabel('y')
title('Regione di assoluta stabilità R-K 1 stadio')

Z1 = min(Z, 1);

figure(2)
mesh(X, Y, Z1)
xlabel('x')
ylabel('y')
zlabel('z')
title('Regione di assoluta stabilità R-K 1 stadio')

%% Secondo metodo

[X, Y] = meshgrid(x);

Z = abs(1-X-Y.*1i);

v = [1 1];

figure(3)
contourf(X, Y, Z, v)
grid on
grid minor
xlabel('x')
ylabel('y')
title('Regione di assoluta stabilità R-K 1 stadio implicito')

Z1 = max(Z, 1);

figure(4)
mesh(X, Y, Z1)
xlabel('x')
ylabel('y')
zlabel('z')
title('Regione di assoluta stabilità R-K 1 stadio implicito')

%% Terzo metodo

s = 2;

alpha = 1/2;
% alpha = 1;

A = [0 0; alpha 0];
b = [1-1/2/alpha 1/2/alpha];

[X, Y] = meshgrid(x);

I = eye(s);

r =@(x,y) abs(1+(x+y.*1i)*b*((I-(x+y.*1i)*A)\ones(s,1)));

F=zeros(n);
for i=1:n
    for j=1:n
        F(i,j) = r(X(i,j),Y(i,j));
    end
end

v = [1 1];

figure(5)
contour(X, Y, F, v)
grid on
grid minor
xlabel('x')
ylabel('y')
title('Regione di assoluta stabilità R-K 2 stadi')

F1 = min(F, 1);

figure(6)
mesh(X, Y, F1)
xlabel('x')
ylabel('y')
zlabel('z')
title('Regione di assoluta stabilità R-K 2 stadi')

%% Quarto metodo

n = 100;
x = linspace(-4,4,n);
y = x;

A = [0 0 0 0; 1/2 0 0 0 ; 0 1/2 0 0; 0 0 1 0];
b = [1/6 1/3 1/3 1/6];

s = 4;

[X, Y] = meshgrid(x);

I = eye(s);

r =@(x,y) abs(1+(x+y.*1i)*b*((I-(x+y.*1i)*A)\ones(s,1)));

F=zeros(n);
for i=1:n
    for j=1:n
        F(i,j) = r(X(i,j),Y(i,j));
    end
end

v = [1 1];

figure(7)
contour(X, Y, F, v)
grid on
grid minor
xlabel('x')
ylabel('y')
title('Regione di assoluta stabilità R-K 4 stadi')

F1 = min(F, 1);

figure(8)
mesh(X, Y, F1)
xlabel('x')
ylabel('y')
zlabel('z')
title('Regione di assoluta stabilità R-K 4 stadi')








