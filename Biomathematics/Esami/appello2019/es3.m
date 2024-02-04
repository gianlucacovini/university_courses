% Fatto

clear
close all

%% Dati

% Iapp = @(t) 2*(t<=0.1);
Iapp = @(t) 0;
b = 5;
c = 1;
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

tspan = [0 200];

%% Solver

dF_RMC = @(t, x) [b*x(1)*(x(1)-beta)*(delta-x(1))-c*x(1)*x(2)+Iapp(t);...
                e*(x(1)-gamma*x(2))];

dF_FHN = @(t, x) [b*x(1)*(x(1)-beta)*(delta-x(1))-c*x(2)+Iapp(t);...
                e*(x(1)-gamma*x(2))];

v0 = 0.3;
w0 = 0.1;


% v0 = -0.1;
% w0 = -0.1;

options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[t, x_RMC] = ode15s(dF_RMC, tspan, [v0; w0], options);
[~, x_FHN] = ode15s(dF_FHN, t, [v0; w0], options);

v_RMC = x_RMC(:, 1);
w_RMC = x_RMC(:, 2);

v_FHN = x_FHN(:, 1);
w_FHN = x_FHN(:, 2);

%% Direction field

v_values_df = linspace(-0.3, 1.3, 20); % Adjust as needed
w_values_df = linspace(-0.2, 1.3, 20); % Adjust as needed
[V_df, W_df] = meshgrid(v_values_df, w_values_df);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_df = b*V_df.*(V_df - beta).*(delta - V_df) - c*V_df.*W_df + Iapp(t);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_df = e*(V_df - gamma*W_df);  % Change in w

%% Nullclines FHN

v_values_nc = linspace(-1, 2, 1000); % Adjust as needed
w_values_nc = linspace(-1, 2, 1000); % Adjust as needed
[V_nc, W_nc] = meshgrid(v_values_nc, w_values_nc);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_ncFHN = b*V_nc.*(V_nc - beta).*(delta - V_nc) - c*W_nc + Iapp(t);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_ncFHN = e*(V_nc - gamma*W_nc);  % Change in w

%% Nullclines RMC

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_ncRMC = b*V_nc.*(V_nc - beta).*(delta - V_nc) - c*V_nc.*W_nc + Iapp(t);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_ncRMC = e*(V_nc - gamma*W_nc);  % Change in w

%% Plots

figure(1)
hold on
contour(V_nc, W_nc, DV_ncFHN, [0 0], '--', 'LineWidth', 2.5)
contour(V_nc, W_nc, DW_ncFHN, [0 0], '--', 'LineWidth', 2.5)
contour(V_nc, W_nc, DV_ncRMC, [0 0])
contour(V_nc, W_nc, DW_ncRMC, [0 0])
hold off
% Controllare la linea verticale che non mi convince

figure(2)
subplot(1, 2, 1)
hold on
plot(t, v_RMC, t, w_RMC)
plot(t, v_FHN, '--')
plot(t, w_FHN, '--')
hold off
grid on
ylim([-0.5 1.2])
subplot(1, 2, 2)
hold on
plot(v_RMC, w_RMC, 'LineWidth', 1.5)
plot(v_FHN, w_FHN, 'k--', 'LineWidth', 1.5)
quiver(V_df, W_df, DV_df, DW_df);
contour(V_nc, W_nc, DV_ncFHN, [0 0], 'r--', 'LineWidth', 1.5)
contour(V_nc, W_nc, DW_ncFHN, [0 0], 'r--', 'LineWidth', 1.5)
contour(V_nc, W_nc, DV_ncRMC, [0 0])
contour(V_nc, W_nc, DW_ncRMC, [0 0])
hold off
axis([-0.3 1.3 -0.2 1.3])
xlabel('v')
ylabel('w')