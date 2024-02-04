clear
close all
clc

I=[0 1];
IT=[0 10];

n=100; m=300;

alpha=.1;
sig=1e-3;
b=5;
u0=@(x) zeros(size(x));



f=@(u) b.*u.*(1-u); % dato di KF
u_x0=0; u_x1=0; %dato al bordo di Neumann

[U,X,T]=KF_solver(f,I,IT,n,m,u_x0,u_x1,u0,sig,alpha);
%Iapp chiamata a parte quando serve

[T_mesh, X_mesh] = meshgrid(T, X);

% proviamo il plot animato

% figure(1)
% for i=2:length(T)
%         mesh(T(1:i), X, U(:,1:i))
%         axis([0 10  0 1  0 1])
%         drawnow
% %         pause(0.1)
% end

figure(2)
subplot(3, 1, 1)
plot(T, U(floor(n/2), :))
ylim([-0.5 1.5])
xlabel('t')
ylabel('x')
title('v(x=0.5, t)')

subplot(3, 1, 2)
imagesc(X, T, U')
set(gca, 'YDir', 'normal');
colorbar
xlabel('x')
ylabel('t')
title('imagesc v(x, t)')

subplot(3, 1, 3)
surf(X, T, U')
colorbar
xlabel('x')
ylabel('t')
title('surf v(x, t)')


for i=2:length(T)
        figure(3)
        plot(X, U(:,i), 'r-o')
        axis([0 1  -0.2 1.5])
        title('t =', T(i))
        xlabel('x')
        ylabel('u')
        drawnow
%         pause(0.1)
end

