clc
clear all
% close all

%% Parameters

Cm = 1; %[muF/cm2]
gNa = 120; %[mOhm-1/cm2]
gK = 36; %[mOhm-1/cm2]
gL = 0.3; %[mOhm-1/cm2]
VNa = 115; %[mV]
VK = -12; %[mV]
VL = 10.6; %[mV]

T = 50;

%% Dati iniziali di resting
% V0 = 2.7570e-4; %-0.195678822932483;%2.7570e-4;
% n0 = 3.1768e-1;
V0 = 2.7570e-4;
n0 = 3.1768e-01; 
%se do i valori di resting mi accorgo dall'orbita che non parte il
%potenziale.

%% Iapp
Iapp = 10;

Tstim = 1; %se Tstim==T aggiorno i plot

%% Solver

opts = [];
[tt,yy] = ode15s(@HHkinetik,[0 T],[V0;n0],opts,Cm,gNa,gK,gL,VNa,VK,VL,Iapp,Tstim);
v = yy(:,1);
n = yy(:,2);

%

figure()
comet(tt,v)
grid on

% Nullclines
[V,N] = meshgrid(-20:0.1:120, 0:0.01:1); %definisco i punti su cui voglio vedere le nullc  
%se non faccio la grid avrei 2 vettori di dimensioni diverse che poi si
%moltiplicano in nullc_V perciò non verrà mai!!!

alfa_m = 0.1*(25-V).*(exp((25-V)/10)-1).^(-1);
alfa_n = 0.01*(10-V).*(exp((10-V)/10)-1).^-1;
beta_m = 4*exp(-V/18);
beta_n = 0.125*exp(-V/80);

tau_n = 1./(alfa_n+beta_n);
tau_m = 1./(alfa_m+beta_m);
n_inf = alfa_n./(alfa_n+beta_n);
m_inf = alfa_m./(alfa_m+beta_m);

if Tstim == T
    nullc_V2 = -(gNa*m_inf.^3.*(0.8-N).*(V-VNa) + gK*N.^4.*(V-VK) + gL*(V-VL)) + Iapp;
end

    nullc_V = -(gNa*m_inf.^3.*(0.8-N).*(V-VNa) + gK*N.^4.*(V-VK) + gL*(V-VL));
    nullc_n = n_inf;

figure()
contour(V,N,nullc_V, [0 0],'r')
hold on
%     contour(V,n,nullc_V2, [0 0])
%     hold on
plot(V(1,:),nullc_n(1,:))
xlabel('V')
ylabel('n')
comet(v,n)
legend('nullc-V','nullc-m','traiett')%,'nullc-V-Iapp'


%% ES. 3 ESAME 15/01/21

%calcolare punto P, minimo del pot d'azione

v_p=min(v);
t_p=tt(find(v==v_p));
v_picco1 = max(v);

t_o=15; %Da grafico

%% ES. 3 ESAME 15/01/21

%secondo stimolo tra t_p e t_o
Iapp=10;
Tstim=1;

Iapp2=38;
Tstim2=Tstim; %sempre

%t_bar è il tempo compreso tra t_p e t_o
t_bar=[5,9,13]; %prendo come valori 5,9,13(compresi tra 3,2 e 20)

opts = [];

for k = 1:length(t_bar)
[tt,yy] = ode15s(@HHkinetik2,[0 T],[V0;n0],opts,Cm,gNa,gK,gL,VNa,VK,VL,Iapp,Iapp2,Tstim,Tstim2,t_bar(k));
v = yy(:,1);
n = yy(:,2);

figure(3)
plot(tt,v)
grid on
title('Potenziale nel tempo')
hold on

[V,N] = meshgrid(-20:0.1:120, 0:0.01:1);

alfa_m = 0.1*(25-V).*(exp((25-V)/10)-1).^(-1);
alfa_n = 0.01*(10-V).*(exp((10-V)/10)-1).^-1;
beta_m = 4*exp(-V/18);
beta_n = 0.125*exp(-V/80);

tau_n = 1./(alfa_n+beta_n);
tau_m = 1./(alfa_m+beta_m);
n_inf = alfa_n./(alfa_n+beta_n);
m_inf = alfa_m./(alfa_m+beta_m);

nullc_V = -(gNa*m_inf.^3.*(0.8-N).*(V-VNa) + gK*N.^4.*(V-VK) + gL*(V-VL));
nullc_n = n_inf;
% nullc_V = -(gK*N.^4.*VK+gNa*m_inf.^3.*(0.8-N).*VNa+gL*VL)/(gK*N.^4+gNa*m_inf.^3.*(0.8-N)+gL);

figure(4)
contour(V,N,nullc_V, [0 0],'r')
hold on
plot(V(1,:),nullc_n(1,:))
xlabel('V')
ylabel('N')
comet(v,n)
legend('nullc-V','nullc-m','traiett')
hold on

pause
end