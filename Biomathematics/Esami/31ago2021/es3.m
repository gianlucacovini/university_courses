clear
close all
clc

%% FHN1D

I=[0 1];
IT=[0 10];

h=0.01; delta_t=0.05;

sigma = 0.01; 
b = 4; 
c = 1; 
beta = 0.1; 
delta = 1; 
gamma = 0.25; 
e = 0.1;


v0=@(x) zeros(size(x));
w0=@(x) zeros(size(x));

f=@(v, w) b.*v.*(v-beta).*(delta-v)-c*w; % dato di FHN1D

% [V, W, X, T]=FHN1D_solver2_periodic(f,I,IT,h,delta_t,v0,w0,sigma,gamma,e);
[V, W, X, T]=FHN1D_solver_periodic_puntob(f,I,IT,h,delta_t,v0,w0,sigma,gamma,e);


V = V';
W = W';

n = length(X);
m = length(T);

figure(2)
subplot(3, 1, 1)
plot(T, V(floor(n/2), :))
ylim([-0.5 1.5])
xlabel('t')
ylabel('x')
title('v(x=0.5, t)')

subplot(3, 1, 2)
imagesc(X, T, V')
set(gca, 'YDir', 'normal');
colorbar
xlabel('x')
ylabel('t')
title('imagesc v(x, t)')

subplot(3, 1, 3)
surf(X, T, V')
colorbar
xlabel('x')
ylabel('t')
title('surf v(x, t)')


for i=2:length(T)
        figure(3)
        plot(X, V(:,i), 'r-o')
        axis([0 1 -0.3 1.3])
        title('t =', T(i))
        xlabel('x')
        ylabel('u')
        drawnow
%         pause(0.1)
end

