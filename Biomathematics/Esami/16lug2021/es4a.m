clear
close all
clc

I=[0 1];
IT=[0 50];

Sig = 0.005;
C_m = 1;
G_na = 120;
G_k = 36;
G_l = 0.3;
V_na = 115;
V_k = -12;
V_l = 10.6;

alpha_m=@(v) 0.1.*(25-v)./(exp((25-v)./10)-1);
alpha_h=@(v) 0.07.*exp(-v./20);
alpha_n=@(v) 0.01.*(10-v)./(exp((10-v)./10)-1);
beta_m=@(v) 4.*exp(-v./18);
beta_h=@(v) 1./(exp((30-v)./10)+1);
beta_n=@(v) 0.125.*exp(-v./80);

f=@(V, m, h, n) -G_na.*(m.^3).*h.*(V-V_na)-G_k.*(n.^4).*(V-V_k)-G_l.*(V-V_l);

% ATTENZIONE: i valori di resting non sono 0
v0=@(x) 2.7570e-04*ones(size(x));
m0=@(x) 5.2934e-02*ones(size(x));
h0=@(x) 5.9611e-01*ones(size(x));
n0=@(x) 3.1768e-01*ones(size(x));

v_x0=0; v_x1=0; %dato al bordo di Neumann

% [V, m, h, n, X,T]=HH1D_solver(f,I,IT,n_x,n_t,v_x0,v_x1,v0,m0,h0,n0,Sig,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n);
%Iapp chiamata a parte quando serve
[V, m, h, n, X,T]=HH1D_solver2_periodic(f,I,IT,0.01,0.05,v_x0,v_x1,v0,m0,h0,n0,Sig,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n);

X = X';
T = T';
V = V';
m = m';
n = n';
h = h';

g_na=G_na.*(m.^3).*h;
g_k=G_k.*n.^4;

I_k=g_k.*(V-V_k);
I_na=g_na.*(V-V_na);
I_l=G_l.*(V-V_l);
I_ion=I_na+I_k+I_l;

[n_x, ~] = size(X);

% figure(1)
% for i=2:length(T)
%         subplot(3, 1, 1)
%         plot(X, V(:,i), 'o-')
%         axis([0 1  -20 120])
%         drawnow
%         
%         subplot(3, 1, 2)
%         plot(X, m(:,i), X, h(:,i), X, n(:,i))
%         axis([0 1  -0.5 1])
%         legend('m', 'h', 'n')
% 
%         subplot(3, 1, 3)
%         plot(X, I_k(:, i), X, I_na(:, i), X, I_ion(:, i))
%         axis([0 1 -1000 2500])
%         legend('I_K', 'I_Na', 'I_ion')
% %         pause(0.1)
% end

figure(1)
subplot(3, 2, 1)
plot(T, V(20, :))
title('v')

subplot(3, 2, 2)
plot(T, m(20, :))
title('m')

subplot(3, 2, 3)
plot(T, h(20, :))
title('h')

subplot(3, 2, 4)
plot(T, n(20, :))
title('n')

subplot(3, 2, 5)
plot(T, I_k(20, :), T, I_na(20, :))
legend('I_K', 'I_Na')

subplot(3, 2, 6)
plot(T, I_ion(20, :))
title('Itot')


index = find(abs(T-5.05)<0.01);

figure(2)
subplot(3, 2, 1)
plot(X, V(:, index))
title('v')

subplot(3, 2, 2)
plot(X, m(:, index))
title('m')

subplot(3, 2, 3)
plot(X, h(:, index))
title('h')

subplot(3, 2, 4)
plot(X, n(:, index))
title('n')

subplot(3, 2, 5)
plot(X, I_k(:, index), X, I_na(:, index))
legend('I_K', 'I_Na')

subplot(3, 2, 6)
plot(X, I_ion(:, index))
title('Itot')