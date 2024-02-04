clear
close all
clc

I=[0 0.5];
IT=[0 25];

Sig = 0.001;
C_m = 1;
% G_na = @(x, t) 120*(x>0.05) + 120*(t>4).*(x<=0.05);
G_na = @(x, t) 120;
% G_k = @(x, t) 36;
G_k = @(x, t) 36*(x>0.05) + 36*(t>4).*(x<=0.05);
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

f=@(V, m, h, n, x, t) -G_na(x, t)'.*(m.^3).*h.*(V-V_na)-G_k(x, t)'.*(n.^4).*(V-V_k)-G_l.*(V-V_l);

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

for k=1:length(T)
    G_na_matrix(:, k) = G_na(X, T(k));
end

for k=1:length(T)
    G_k_matrix(:, k) = G_k(X, T(k));
end

g_na=G_na_matrix.*(m.^3).*h;
g_k=G_k_matrix.*n.^4;

I_k=g_k.*(V-V_k);
I_na=g_na.*(V-V_na);
I_l=G_l.*(V-V_l);
I_ion=I_na+I_k+I_l;

[n_x, ~] = size(X);

figure(1)
for i=2:length(T)
        plot(X, V(:,i), 'o-')
        axis([0 0.5  -20 120])
        drawnow
        
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
end

index = find(X == 0.25);

figure(2)
subplot(3, 2, 1)
plot(T, V(index, :))
xlabel('t')
ylabel('V')
title('v')

subplot(3, 2, 2)
plot(T, m(index, :))
xlabel('t')
ylabel('m')
title('m')

subplot(3, 2, 3)
plot(T, h(index, :))
xlabel('t')
ylabel('h')
title('h')

subplot(3, 2, 4)
plot(T, n(index, :))
xlabel('t')
ylabel('n')
title('n')

subplot(3, 2, 5)
plot(T, I_k(index, :), T, I_na(index, :))
xlabel('t')
ylabel('I')
legend('I_K', 'I_{Na}', 'Location', 'bestoutside')
title('I_k e I_ {Na}')

subplot(3, 2, 6)
plot(T, I_ion(index, :))
xlabel('t')
ylabel('I_ion')
title('Itot')

sgtitle('Dinamica del punto x=0.25 - Potassio bloccato')


index = find(T == 3);

figure(3)
subplot(3, 2, 1)
plot(X, V(:, index))
xlabel('t')
ylabel('V')
title('v')

subplot(3, 2, 2)
plot(X, m(:, index))
xlabel('t')
ylabel('m')
title('m')

subplot(3, 2, 3)
plot(X, h(:, index))
xlabel('t')
ylabel('h')
title('h')

subplot(3, 2, 4)
plot(X, n(:, index))
xlabel('t')
ylabel('n')
title('n')

subplot(3, 2, 5)
plot(X, I_k(:, index), X, I_na(:, index))
xlabel('t')
ylabel('I')
legend('I_K', 'I_{Na}', 'Location', 'bestoutside')
title('I_K e I_{Na}')

subplot(3, 2, 6)
plot(X, I_ion(:, index))
xlabel('t')
ylabel('I_ion')
title('Itot')

sgtitle('Dinamica al tempo t=3 - Potassio bloccato')

index = find(T == 6);

figure(4)
subplot(3, 2, 1)
plot(X, V(:, index))
xlabel('t')
ylabel('V')
title('v')

subplot(3, 2, 2)
plot(X, m(:, index))
xlabel('t')
ylabel('m')
title('m')

subplot(3, 2, 3)
plot(X, h(:, index))
xlabel('t')
ylabel('h')
title('h')

subplot(3, 2, 4)
plot(X, n(:, index))
xlabel('t')
ylabel('n')
title('n')

subplot(3, 2, 5)
plot(X, I_k(:, index), X, I_na(:, index))
xlabel('t')
ylabel('I')
legend('I_K', 'I_{Na}', 'Location', 'bestoutside')
title('I_K e I_{Na}')

subplot(3, 2, 6)
plot(X, I_ion(:, index))
xlabel('t')
ylabel('I_ion')
title('Itot')

sgtitle('Dinamica al tempo t=6 - Potassio bloccato')