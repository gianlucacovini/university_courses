% Tutto ok tranne I_ion!! :(

clear
close all

I=[0 1];
IT=[0 40];

values = [40, 60, 80, 100, 120, 140];

count = 1;

for G_na=values

Sig = 1e-3;
C_m = 1;
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
[V, m, h, n, X,T]=HH1D_solver2(f,I,IT,0.01,0.05,v_x0,v_x1,v0,m0,h0,n0,Sig,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n);

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
%         plot(X, V(:,i))
%         axis([0 1  -20 120])
%         title("t = ", T(i))
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

%% Velocità

[~, indice_in] = max(V(1, :));
[~, indice_fin] = max(V(end, :));
t_prop(count) = T(indice_fin)-T(indice_in);
vel_potenziale(count) = 1/t_prop(count); %[cm/msec]

count = count +1;

end

figure(2)
plot(values, vel_potenziale, '-o')


