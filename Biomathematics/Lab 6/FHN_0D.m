% Fare animazioni

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

tspan = [0 100];

%% Solver

dF = @(t, x) [b*x(1)*(x(1)-beta)*(delta-x(1))-c*x(2)+Iapp(t);...
                e*(x(1)-gamma*x(2))];

v0 = 0.6;
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
DV_df = b*V_df.*(V_df - beta).*(delta - V_df) - c*W_df + Iapp(t);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_df = e*(V_df - gamma*W_df);  % Change in w

%% Nullclines

v_values_nc = linspace(-1, 2, 1000); % Adjust as needed
w_values_nc = linspace(-1, 2, 1000); % Adjust as needed
[V_nc, W_nc] = meshgrid(v_values_nc, w_values_nc);

% Calculate the direction (dv/dt, dw/dt) at each grid point
DV_nc = b*V_nc.*(V_nc - beta).*(delta - V_nc) - c*W_nc + Iapp(t);  % Mettere Iapp(0.5) ad esempio nel caso di impulso iniziale
DW_nc = e*(V_nc - gamma*W_nc);  % Change in w

%% Plots

figure(1)
subplot(1, 2, 1)
plot(t, v, t, w)
grid on
ylim([-0.5 1.2])
subplot(1, 2, 2)
hold on
plot(v, w)
quiver(V_df, W_df, DV_df, DW_df);
contour(V_nc, W_nc, DV_nc, [0 0])
contour(V_nc, W_nc, DW_nc, [0 0])
hold off
axis([-0.4 1 -0.1 0.8])
xlabel('v')
ylabel('w')

% for i=2:length(T)
%         figure(3)
%         plot(X, U(:,i), 'r-o')
%         axis([0 1  -0.2 1.5])
%         title('t =', T(i))
%         xlabel('x')
%         ylabel('u')
%         drawnow
% %         pause(0.1)
% end

%% OSSERVAZIONI

% Si vede l'effetto soglia. Il potenziale d'azione per v0 > di una certa
% quantità va ad un picco

% Si vedono le frecce. Abbiamo per v0=0.6 che il punto tende al punto di
% equilibrio 0,0 che è intersezione delle nullclines. è equilibrio stabile
% perché sta nel ramo h- (visto a lezione). Prima di arrivare al pto di
% equilibrio fa un lungo excursus perché siamo sopra soglia. Va dalla prima
% nullcline e lì il campo è solo verticale, poi quando interseca l'altra è
% solo orizzontale.
% Per v0=0.1 siamo sotto soglia e niente grande excursus. Se sta a sinistra
% della cubica è sotto soglia, niente excursus
% Si può fare la stessa cosa con w. Provare nel caso con w0 = -0.1

% Se aumentiamo Iapp, tipo a 0.5, si genera potenziale d'azione perché la
% cubica si è spostata in alto. Notiamo che il pto di equilibrio
% (intersezione nullclines) cade ora nel ramo instabile (h0), quindi non è
% più punto di equilibrio ma abbiamo ciclo limite.

% Aumentando Iapp va sempre più su la cubica e segue la stessa dinamica,
% fino a 2.5 dove abbiamo che il pto di equilibrio sta in h+ che è ancora
% stabile

% Con l'impulso invece, se è basso non c'è potenziale d'azione (rimane
% sotto soglia). Se invece è più alto parte il potenziale d'azione.
% L'impulso deve portare oltre la zona soglia (oltre la cubica). Si può
% aumentare la durata o l'intensità dello stimolo nel caso



% Attenzione: le nullcline per un impulso di corrente sono quelle dopo
% l'impulso (quindi per I = 0)










