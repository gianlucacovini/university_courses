clear
close all
clc

%% Dati

Iapp = @(t) 5*(t<=0.1);
b = 1;
c = 1;
beta = 0.1;
delta = 1;
gamma = 0.5;
e = 0.01;

tspan = [0 150];

%% Solver

dF = @(t, x) [b*x(1)*(x(1)-beta)*(delta-x(1))-c*x(2)+Iapp(t);...
                e*(x(1)-gamma*x(2))];

v0 = 0;
w0 = 0;

options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[t, x] = ode15s(dF, tspan, [v0; w0], options);

v = x(:, 1);
w = x(:, 2);

[n, ~] = size(v);

%% Direction field

v_values_df = linspace(-0.4, 1, 20); % Adjust as needed
w_values_df = linspace(-0.1, 0.8, 20); % Adjust as needed
[V_df, W_df] = meshgrid(v_values_df, w_values_df);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_df = b*V_df.*(V_df - beta).*(delta - V_df) - c*W_df + Iapp(0.5);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_df = e*(V_df - gamma*W_df);  % Change in w

%% Nullclines

v_values_nc = linspace(-1, 2, 1000); % Adjust as needed
w_values_nc = linspace(-1, 2, 1000); % Adjust as needed
[V_nc, W_nc] = meshgrid(v_values_nc, w_values_nc);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_nc = b*V_nc.*(V_nc - beta).*(delta - V_nc) - c*W_nc + Iapp(0.5);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_nc = e*(V_nc - gamma*W_nc);  % Change in w

%% Plots

figure(1)
subplot(1, 2, 1)
plot(t, v, t, w)
ylim([-0.5 1.5])
xlabel('t')
grid on
subplot(1, 2, 2)
hold on
plot(v, w)
axis([-0.4 1 -0.1 0.8])
quiver(V_df, W_df, DV_df, DW_df);
contour(V_nc, W_nc, DV_nc, [0 0])
contour(V_nc, W_nc, DW_nc, [0 0])
hold off
xlabel('v')
ylabel('w')

%% Durata potenziale d'azione (DA VERIFICARE)

% è il tempo per tornare al potenziale di resting a quanto ho capito
[~, index_min] = min(v);
tmin_v = t(index_min); 

t0 = t(abs(v-v0)<1e-3);

indice0 = min(find(t0>tmin_v));
t_fin = t0(indice0);
[~, index_max] = max(v); 
t_in = t(index_max);

APD = t_fin - t_in