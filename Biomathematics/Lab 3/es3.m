%% Cooperatività positiva
clear
close all

tspan = [0 100];

k1 = 0.5010;
k_m1 = 500;
k2 = 1;
k3 = 502000;
k_m3 = 500;
k4 = 2;

e0 = 0.1;

s0 = 5;
c1_0 = 0;
c2_0 = 0;
y0 = [s0 c1_0 c2_0];

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t_pos, y_pos] = ode15s(@(t, y) [-k1*y(1)*(e0-y(2)-y(3))+k_m1*y(2)-k3*y(1)*y(2)+k_m3*y(3); ...
                        k1*y(1)*(e0-y(2)-y(3))-(k_m1+k2)*y(2)-k3*y(1)*y(2)+(k4+k_m3)*y(3);...
                        k3*y(1)*y(2)-(k4+k_m3)*y(3)],...
                tspan,...
                y0,...
                options);

s_pos = y_pos(:, 1);
c1_pos = y_pos(:, 2);
c2_pos = y_pos(:, 3);

K1 = (k_m1+k2)/k1;
K2 = (k_m3+k4)/k3;

options_qs = odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);
[~, s_qs_pos] = ode15s(@(t, s) (-e0*s*(k2*K2+k4*s)/(K1*K2+K2*s+s^2)),...
                t_pos,...
                s0,...
                options_qs);

c1_qs_pos = e0*K2*s_qs_pos./(K1*K2+K2*s_qs_pos+s_qs_pos.^2);
c2_qs_pos = e0*s_qs_pos.^2./(K1*K2+K2*s_qs_pos+s_qs_pos.^2);

V_pos = k2*c1_pos + k4*c2_pos;

k_m_squared = (k_m1+k2)*(k_m3+k4)/k1/k3;
V_pos_qs = e0*k4*s_qs_pos.^2./(k_m_squared + s_qs_pos.^2);

%% Cooperatività indipendente

k1 = 102;
k_m1 = 50;
k2 = 1;
k3 = 26;
k_m3 = 50;
k4 = 2;

s0 = 5;
c1_0 = 0;
c2_0 = 0;
y0 = [s0 c1_0 c2_0];

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t_ind, y_ind] = ode15s(@(t, y) [-k1*y(1)*(e0-y(2)-y(3))+k_m1*y(2)-k3*y(1)*y(2)+k_m3*y(3); ...
                        k1*y(1)*(e0-y(2)-y(3))-(k_m1+k2)*y(2)-k3*y(1)*y(2)+(k4+k_m3)*y(3);...
                        k3*y(1)*y(2)-(k4+k_m3)*y(3)],...
                tspan,...
                y0,...
                options);

s_ind = y_ind(:, 1);
c1_ind = y_ind(:, 2);
c2_ind = y_ind(:, 3);

K1 = (k_m1+k2)/k1;
K2 = (k_m3+k4)/k3;

options_qs = odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);
[~, s_qs_ind] = ode15s(@(t, s) (-e0*s*(k2*K2+k4*s)/(K1*K2+K2*s+s^2)),...
                t_ind,...
                s0,...
                options_qs);

c1_qs_ind = e0*K2*s_qs_ind./(K1*K2+K2*s_qs_ind+s_qs_ind.^2);
c2_qs_ind = e0*s_qs_ind.^2./(K1*K2+K2*s_qs_ind+s_qs_ind.^2);

V_ind = k2*c1_ind + k4*c2_ind;

K = 2*K1;
V_ind_qs = 2*k2*e0*s_qs_ind./(K+s_qs_ind);

%% Cooperatività negativa

k1 = 1002;
k_m1 = 500;
k2 = 1;
k3 = 5.02;
k_m3 = 500;
k4 = 2;

s0 = 5;
c1_0 = 0;
c2_0 = 0;
y0 = [s0 c1_0 c2_0];

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t_neg, y_neg] = ode15s(@(t, y) [-k1*y(1)*(e0-y(2)-y(3))+k_m1*y(2)-k3*y(1)*y(2)+k_m3*y(3); ...
                        k1*y(1)*(e0-y(2)-y(3))-(k_m1+k2)*y(2)-k3*y(1)*y(2)+(k4+k_m3)*y(3);...
                        k3*y(1)*y(2)-(k4+k_m3)*y(3)],...
                tspan,...
                y0,...
                options);

s_neg = y_neg(:, 1);
c1_neg = y_neg(:, 2);
c2_neg = y_neg(:, 3);

K1 = (k_m1+k2)/k1;
K2 = (k_m3+k4)/k3;

options_qs = odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);
[~, s_qs_neg] = ode15s(@(t, s) (-e0*s*(k2*K2+k4*s)/(K1*K2+K2*s+s^2)),...
                t_neg,...
                s0,...
                options_qs);

c1_qs_neg = e0*K2*s_qs_neg./(K1*K2+K2*s_qs_neg+s_qs_neg.^2);
c2_qs_neg = e0*s_qs_neg.^2./(K1*K2+K2*s_qs_neg+s_qs_neg.^2);

V_neg = k2*c1_neg + k4*c2_neg;

V_neg_qs = e0*s_qs_neg.*(k2+k4/K2*s_qs_neg)./(K1+s_qs_neg+s_qs_neg.^2/K2);

%% Plot

figure(1)
subplot(3, 1, 1)
plot(t_pos, s_pos, 'r', t_ind, s_ind, 'b', t_neg, s_neg, 'k', t_pos, s_qs_pos, 'r--', t_ind, s_qs_ind, 'b--', t_neg, s_qs_neg, 'k--')
xlabel('t')
ylabel('s')
xlim([0 100])
subplot(3, 1, 2)
plot(t_pos, c1_pos, 'r', t_ind, c1_ind, 'b', t_neg, c1_neg, 'k', t_pos, c1_qs_pos, 'r--', t_ind, c1_qs_ind, 'b--', t_neg, c1_qs_neg, 'k--')
xlabel('t')
ylabel('c1')
xlim([0 100])
subplot(3, 1, 3)
plot(t_pos, c2_pos, 'r', t_ind, c2_ind, 'b', t_neg, c2_neg, 'k', t_pos, c2_qs_pos, 'r--', t_ind, c2_qs_ind, 'b--', t_neg, c2_qs_neg, 'k--')
xlabel('t')
ylabel('c2')
xlim([0 100])

figure(2)
plot(s_pos, V_pos, 'r', s_ind, V_ind, 'b', s_neg, V_neg, 'k',...
    s_qs_pos, V_pos_qs, 'r--', s_qs_ind, V_ind_qs, 'b--', s_qs_neg, V_neg_qs, 'k--')
xlim([-1, 5])
%ylim([0, 2])
legend('pos', 'ind', 'neg', 'qs pos', 'qs ind', 'qs neg')
xlabel('s')
title('velocity')

%% Errori

figure(3)
subplot(3, 1, 1)
plot(t_pos, abs(s_pos-s_qs_pos), 'r', t_ind, abs(s_ind-s_qs_ind), 'b', t_neg, abs(s_neg-s_qs_neg), 'k')
xlabel('t')
ylabel('s')
xlim([0 100])
subplot(3, 1, 2)
plot(t_pos, abs(c1_pos-c1_qs_pos), 'r', t_ind, abs(c1_ind-c1_qs_ind), 'b', t_neg, abs(c1_neg-c1_qs_neg), 'k')
xlabel('t')
ylabel('c1')
xlim([0 100])
subplot(3, 1, 3)
plot(t_pos, abs(c2_pos-c2_qs_pos), 'r', t_ind, abs(c2_ind-c2_qs_ind), 'b', t_neg, abs(c2_neg-c2_qs_neg), 'k')
xlabel('t')
ylabel('c2')
xlim([0 100])

%% OSSERVAZIONI

% Sol. completa ha gli strati limite nella fase transitoria. Penso intenda
% quando esplodono all'inizio (si vede zoomando sui primi 1 ms)

% Gli strati limite si vedono anche nel decadimento limite della velocità
% La velocità positiva qs è molto diversa da quella normale

% Con e0 = 0.1 va più lenta perché c'è meno enzima. L'approssimazione è
% migliore