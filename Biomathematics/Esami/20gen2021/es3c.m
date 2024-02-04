clear
close all
clc

%% Dati

st_times = [20 25 30 35 40 45 50 55 60 80 100 140];

count = 1;

for T_stim=st_times

Iapp = @(t) 5*(t<=0.1) + 5*(t>=T_stim).*(t<=T_stim+1);
b = 1;
c = 1;
beta = 0.1;
delta = 1;
gamma = 0.5;
e = 0.01;

%% Solver

dF = @(t, x) [b*x(1)*(x(1)-beta)*(delta-x(1))-c*x(2)+Iapp(t);...
                e*(x(1)-gamma*x(2))];

v0 = 0;
w0 = 0;

% Attenzione alle options se aggiungiamo un impulso a 140, dobbiamo allora
% risolvere due problemi: uno fino a 140 e uno dopo
options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[t1, x1] = ode15s(dF, [0 T_stim], [v0; w0], options);

v1 = x1(:, 1);
w1 = x1(:, 2);

[t2, x2] = ode15s(dF, [T_stim 300], [v0; w0], options);

v2 = x2(:, 1);
w2 = x2(:, 2);

v = [v1; v2];
w = [w1; w2];

t = [t1; t2];

%% Direction field

v_values_df = linspace(-0.4, 2.5, 20); % Adjust as needed
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
ylim([-0.5 2.5])
xlabel('t')
grid on
subplot(1, 2, 2)
hold on
plot(v, w)
axis([-0.4 2.5 -0.1 0.5])
quiver(V_df, W_df, DV_df, DW_df);
contour(V_nc, W_nc, DV_nc, [0 0])
contour(V_nc, W_nc, DW_nc, [0 0])
hold off
xlabel('v')
ylabel('w')

%% Durata potenziale d'azione (DA VERIFICARE)

[~, index_max] = max(v); 
t_in = t(index_max);

[~, index_min] = min(v(index_max:end));
tmin_v = t(index_max+index_min); 

t0 = t(abs(v-v0)<1e-3);

indice0 = min(find(t0>tmin_v));
t_fin = t0(indice0);

APD(count) = t_fin % - t_in % Con o senza - t_in ???

count = count + 1;

end

%% APD plot

figure(2)
plot(st_times, APD)
% ylim([0 150])

