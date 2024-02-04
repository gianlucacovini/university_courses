clear
close all
clc

% Consideriamo il modello completo e poi proiettiamo le soluzioni sullo
% spazio veloce-lento

%% dati
C_m = 1;
G_na = 120;
G_k = 36;
G_l = 0.3;
V_na = 115;
V_k = -12;
V_l = 10.6;

for alpha=-19:-0.1:-20

Iapp = @(t) alpha*(t <= 1);

a_m=@(v) 0.1.*(25-v)./(exp((25-v)./10)-1);
a_h=@(v) 0.07.*exp(-v./20);
a_n=@(v) 0.01.*(10-v)./(exp((10-v)./10)-1);
b_m=@(v) 4.*exp(-v./18);
b_h=@(v) 1./(exp((30-v)./10)+1);
b_n=@(v) 0.125.*exp(-v./80);

I=[0 50];

%% variabili V,m,h,n --> fai function

F=@(t,x) [-G_na.*x(3).*(x(2).^3).*(x(1)-V_na)-G_k.*(x(4).^4).*(x(1)-V_k)-G_l.*(x(1)-V_l)+Iapp(t);...
            a_m(x(1)).*(1-x(2))-b_m(x(1))*x(2);...
            a_h(x(1)).*(1-x(3))-b_h(x(1))*x(3);...
            a_n(x(1)).*(1-x(4))-b_n(x(1))*x(4)];

v0=2.7570e-04;
m0=5.2934e-2;
h0=5.9611e-1;
n0=3.1768e-1;


options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[T,X] = ode15s(F,I,[v0;m0;h0;n0],options);

%% altri valori per i plot
V=X(:,1);
m=X(:,2);
h=X(:,3);
n=X(:,4);

g_na=G_na.*(m.^3).*h;
g_k=G_k.*n.^4;

I_k=g_k.*(V-V_k);
I_na=g_na.*(V-V_na);
I_l=G_l.*(V-V_l);
I_ion=I_na+I_k+I_l;

Tau_n= 1./(a_n(V)+b_n(V));
Tau_m= 1./(a_m(V)+b_m(V));
Tau_h= 1./(a_h(V)+b_h(V));

n_inf= a_n(V)./(a_n(V)+b_n(V));
m_inf= a_m(V)./(a_m(V)+b_m(V));
h_inf= a_h(V)./(a_h(V)+b_h(V));

%% Nullclines

v_values_nc = linspace(-100, 150, 200); % Adjust as needed
n_values_nc = linspace(-0.2, 1.2, 200); % Adjust as needed
[V_nc, N_nc] = meshgrid(v_values_nc, n_values_nc);

% Calculate the direction (dv/dt, dn/dt) at each grid point
DV_nc = -G_na.*(0.8-N_nc).*(a_m(V_nc)./(a_m(V_nc)+b_m(V_nc)).^3).*(V_nc-V_na)-G_k.*(N_nc.^4).*(V_nc-V_k)-G_l.*(V_nc-V_l)+Iapp(10);  % Attenzione a mettere Iapp(0.5), deve essere effettivamente 0
DN_nc = a_n(V_nc).*(1-N_nc)-b_n(V_nc).*N_nc;  % Change in w

%% plot

figure(1)
hold on
txt = ['I_{app} =',num2str(alpha)];
plot(T, V, "DisplayName",txt)    
grid on
% ylim([-0.5 1.2])
legend show

figure(2)
hold on
plot(V, n, "DisplayName",txt)
% quiver(V_df, N_df, DV_df, DN_df);
contour(V_nc, N_nc, DV_nc, [0 0])
contour(V_nc, N_nc, DN_nc, [0 0])
axis([-80 130 0 0.8])
xlabel('V')
ylabel('n')

end
hold off
legend show

% La soglia è -19.4
