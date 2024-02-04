clear
close all

I=[0 1];
IT=[0 45];

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

% ATTENZIONE: i valori di resting non sono 0
v0=@(x) 2.7570e-04*ones(size(x));
m0=@(x) 5.2934e-02*ones(size(x));
h0=@(x) 5.9611e-01*ones(size(x));
n0=@(x) 3.1768e-01*ones(size(x));

v_x0=0; v_x1=0; %dato al bordo di Neumann

gNaSalto = 27.625;
% 25.6

%Iapp chiamata a parte quando serve
[V, m, h, n, X,T]=HH1D_solver_es4b(gNaSalto, I,IT,0.01,0.05,v_x0,v_x1,v0,m0,h0,n0,Sig,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n);


X = X';
T = T';
V = V';
m = m';
n = n';
h = h';

figure(1)
for i=2:length(T)
        plot(X, V(:,i))
        axis([0 1  -20 120])
        title("t = ", T(i))
        drawnow
end

