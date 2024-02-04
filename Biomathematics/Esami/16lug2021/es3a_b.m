% Fare animazioni

clear
close all

%% Dati

Iapp = @(t) -100*(t<=0.1);
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
DV_df = b*V_df.*(V_df - beta).*(delta - V_df) - c*W_df + Iapp(10);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_df = e*(V_df - gamma*W_df);  % Change in w

%% Nullclines

v_values_nc = linspace(-3, 3, 1000); % Adjust as needed
w_values_nc = linspace(-0.5, 0.5, 1000); % Adjust as needed
[V_nc, W_nc] = meshgrid(v_values_nc, w_values_nc);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_nc = b*V_nc.*(V_nc - beta).*(delta - V_nc) - c*W_nc + Iapp(10);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_nc = e*(V_nc - gamma*W_nc);  % Change in w

%% Plots

figure(1)
subplot(1, 2, 1)
plot(t, v, t, w)
grid on
% ylim([-0.5 1.2])
subplot(1, 2, 2)
hold on
plot(v, w, '--', 'LineWidth', 1.5)
% quiver(V_df, W_df, DV_df, DW_df);
contour(V_nc, W_nc, DV_nc, [0 0])
contour(V_nc, W_nc, DW_nc, [0 0])
hold off
% axis([-0.4 1 -0.1 0.8])
xlabel('v')
ylabel('w')
