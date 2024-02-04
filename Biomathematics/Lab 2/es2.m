%% Punto 1
clear
close all

tspan = [0 100];

k1 = 102;
k_m1 = 50;
k2 = 1;
k3 = 26;
k_m3 = 50;

e0 = 1;

s0 = 5;
i0 = 1;
c1_0 = 0;
c2_0 = 0;
y0 = [s0 i0 c1_0 c2_0];

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t, y] = ode15s(@(t, y) [-k1*y(1)*(e0-y(3)-y(4))+k_m1*y(3); ...
                        -k3*y(2)*(e0-y(3)-y(4))+k_m3*y(4);...
                        k1*y(1)*(e0-y(3)-y(4))-(k_m1+k2)*y(3);...
                        k3*y(2)*(e0-y(3)-y(4))-k_m3*y(4)],...
                tspan,...
                y0,...
                options);

s = y(:, 1);
i = y(:, 2);
c1 = y(:, 3);
c2 = y(:, 4);

k_m = (k_m1+k2)/k1;
k_i = k_m3/k3;

V_max = k2*e0;

[l, ~] = size(t);
i_qs = ones(l, 1);

options_qs = odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);
[t_qs, s_qs] = ode15s(@(t, s) -V_max*s/(s+k_m*(1+1/k_i)),...
                t,...
                s0,...
                options_qs);

c1_qs = e0*k_i*s_qs./(i_qs*k_m + s_qs*k_i + k_i*k_m);
c2_qs = e0*i_qs./(i_qs+k_i*(1+s_qs/k_m));

figure(1)
plot(t, s, 'r', t, i, 'g', t, c1, 'b', t, c2, 'k', t_qs, s_qs, 'r--', t_qs, i_qs, 'g--', t_qs, c1_qs, 'b--', t_qs, c2_qs, 'k--')
legend('s', 'in', 'c1', 'c2', 's_qs', 'in_qs', 'c1_qs', 'c2_qs')
title('competitive inhibition')
xlabel('time')

%% Errori

E_s = abs(s - s_qs);
E_i = abs(i - i_qs);
E_c1 = abs(c1 - c1_qs);
E_c2 = abs(c2 - c2_qs);

figure(2)
plot(t, E_s, 'r', t, E_i, 'g', t, E_c1, 'b', t, E_c2, 'k')
title('error competitive inhibition')
xlabel('time')
legend('E_s', 'E_i', 'E_{c1}', 'E_{c2}')

%% Velocità

V = k2*c1;

V_qs = V_max*s_qs./(s_qs+k_m*(1+i_qs/k_i));

figure(3)
plot(s, V, 'k', s_qs, V_qs, 'k--')
xlabel('substrate')
ylabel('V')
legend('V', 'V_qs')

%% Velocità 3d

figure(4)
plot3(s, i, V, 'k', s_qs, i_qs, V_qs, 'k--')
grid on
xlabel('substrate s')
ylabel('i')
zlabel('V')

%% Errori velocità

E_v = abs(V - V_qs);

figure(5)
plot(t, E_v)
xlabel('time')
title('error velocity')

%% Inibizione allosterica

tspan = [0 100];

s0 = 5;
x0 = 0;
y0 = 0;
z0 = 0;
i0 = 2;

e0 = 1;

u0 = [s0 x0 y0 z0 i0];

k_m1 = 50;
k1 = 102;
k_m3 = 50;
k3 = 26;
k2 = 1;

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t, u] = ode15s(@(t, u) [-k1*(e0-u(2)-u(3)-u(4))*u(1)+k_m1*u(2)+k_m1*u(4)-k1*u(3)*u(1);...
    -k_m1*u(2)+k1*(e0-u(2)-u(3)-u(4))*u(1)-k3*u(2)*u(5)+k_m3*u(4)-k2*u(2);...
    k3*(e0-u(2)-u(3)-u(4))*u(5)-k_m3*u(3)+k_m1*u(4)-k1*u(3)*u(1);...
    -k_m3*u(4)+k3*u(2)*u(5)-k_m1*u(4)+k1*u(3)*u(1);...
    -k3*(e0-u(2)-u(3)-u(4))*u(5)+k_m3*u(3)-k3*u(2)*u(5)+k_m3*u(4)],...
    tspan,...
    u0,...
    options);

s = u(:, 1);
x = u(:, 2);
y = u(:, 3);
z = u(:, 4);
i = u(:, 5);

figure(6)
plot(t, s, 'r', t, x, 'b', t, y, 'k', t, z, 'g', t, i, 'y')
legend('s_{all}', 'c1_{all}', 'c2_{all}', 'c3_{all}', 'i_{all}')
title('allosteric inhibition')
xlabel('time')

ki = k_m3/k3;
ks = k_m1/k1;

V = k2*x; % Attenzione: usare sempre le formule originali! Non quelle approssimate

figure(7)
plot(s,V)
grid on
title("Velocità in funzione di s")
ylabel('V_{all}')
xlabel('s_{all}')
ylim([0 0.7])

figure(8)
plot3(s,i,V)
xlabel("s_{all}")
ylabel("i_{all}")
zlabel("V_{all}")
grid on
title("Velocità in funzione di s e i")

%% OSSERVAZIONI

% Inibizione competitiva:
% Si nota lo strato limite

% Senza inibitore sta più giù perché il sostrato viene consumato più
% velocemente. Infatti anche la velocità è meggiore


% Inibizione allosterica
% Forte strato limite
