%% Punto 1
clear
close all

tspan = [0 100];

k1 = 4.e+6;
k_m1 = 25;
k2 = 15;

K_m = (k_m1+k2)/k1;

s0 = K_m;
c0 = 0;
y0 = [s0 c0];

eps = 1.e-3;

e0 = eps*K_m;

options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t, y] = ode15s(@(t, y) [k_m1*y(2) - k1*y(1)*(e0 - y(2)); ...
                        k1*y(1)*(e0-y(2))-(k_m1+k2)*y(2)],...
                tspan,...
                y0, ...
                options);

s = y(:, 1);
c = y(:, 2);

T_r = 1/k1/s0
T_p = 1/k1/e0

options_u = odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);
[t_u, s_u] = ode15s(@(t, s_u) -k2*e0*s_u/(K_m + s_u),...
                t,...
                s0,...
                options_u);

c_u = e0*s_u./(K_m + s_u)-e0*s_u/(K_m + s0).*exp(-(K_m+s0)*k1*t_u);

figure(1)
plot(t, s/s0, 'r')
hold on
plot(t, c/e0, 'k') 
plot(t_u, s_u/s0, 'ro')
plot(t_u, c_u/e0, 'ko')
hold off
grid on
legend('s', 'c', 's_u', 'c_u')
title(sprintf('normalized solutions for eps=%f', eps))
xlabel('time')

%% Zoom fase transitoria

transitory_time = 10*eps;

% SOLUZIONE INUTILMENTE LUNGA
% [l, ~] = size(t);
% 
% zoom_t = [];
% zoom_s = [];
% zoom_c = [];
% for i=1:l
%     if t(i) <= transitory_time
%         zoom_t = [zoom_t, t(i)];
%         zoom_s = [zoom_s, s(i)];
%         zoom_c = [zoom_c, c(i)];
%     end
% end
% 
% [l_u, ~] = size(t_u);
% 
% zoom_t_u = [];
% zoom_s_u = [];
% zoom_c_u = [];
% for i=1:l_u
%     if t_u(i) <= transitory_time
%         zoom_t_u = [zoom_t_u, t_u(i)];
%         zoom_s_u = [zoom_s_u, s_u(i)];
%         zoom_c_u = [zoom_c_u, c_u(i)];
%     end
% end
% 
% figure(2)
% plot(zoom_t, zoom_s/s0, 'r', zoom_t, zoom_c/e0, 'k', zoom_t_u, zoom_s_u/s0, 'ro', zoom_t_u, zoom_c_u/e0, 'ko')
% grid on
% legend('s', 'c', 's_u', 'c_u')
% title('zoom on [0, 10eps] of normalized solutions')
% xlabel('time')

% SOLUZIONE RAPIDA
figure(2)
plot(t, s/s0, 'r', t, c/e0, 'k', t_u, s_u/s0, 'ro', t_u, c_u/e0, 'ko')
xlim([0, transitory_time])
grid on
legend('s', 'c', 's_u', 'c_u')
title('zoom on [0, 10eps] of normalized solutions for eps=', eps)
xlabel('time')

%% Plot errori

E_s = abs((s-s_u)/s0); 
E_c = abs((c-c_u)/e0);

figure(3)
semilogy(t, E_s, 'r', t, E_c, 'k')
grid on
legend('abs(s-s_u)', 'abs(c-c_u)')
xlabel('time')

%% Spazio delle fasi

figure(4)
plot(s/s0, c/e0, 'b')
hold on
plot(s_u/s0, c_u/e0, 'bo')
xlim([0 1.2])
ylim([0 0.6])
hold off
legend('(s, c)', '(s_u, c_u)')
title('trajectories phase space (s, c)')
xlabel('s')
ylabel('c')

%% OSSERVAZIONI

% Si nota che per valori di eps più grandi (0.1, 1) l'approssimazione
% uniforme perggiora perché l'approssimazione ha un errore di O(eps).
% Con eps più grande le soluzioni sono più veloci, inoltre

% Per tempi più lunghi s e c tendono a 0