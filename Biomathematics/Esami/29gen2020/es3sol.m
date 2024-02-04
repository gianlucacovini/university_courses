% esercizio 3

%% punto a
clc
clear all 
close all

V0 = 2.7570e-04;
m0 = 5.2934e-2;
h0 = 5.9611e-1;
n0 = 3.1768e-1;

% dati iniziali: 
y0 = [V0,m0,h0,n0];
Iapp = -100; 

Cm = 1;
gNa0 = 120;
gK0 = 36;
gL0 = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;

% T = 200;
% Tstim = 1; % durata dell'impulso di corrente
% tspan = [0,Tstim];
% options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
% [t1,y1] = ode15s(@HH,tspan,y0,options, Iapp, Cm, gNa0, gK0 ,gL0, VNa, VK, VL);
% V = y1(:,1);
% n = y1(:,2);
% tspan2 = [t1(end),T];
% y02 = [V(end), n(end)];
% Iapp2 = 0;
% [t2,y2] = ode15s(@HH,tspan2,y02,options, Iapp2, Cm, gNa0, gK0 ,gL0, VNa, VK, VL);
% V = [V; y2(:,1)];
% n = [n; y2(:,2)];
% t = [t1; t2];


T = 50;
Tstim = 1; % durata dell'impulso di corrente
tspan = [0,Tstim];
options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t1,y1] = ode15s(@HH_completo,tspan,y0,options, Iapp, Cm, gNa0, gK0 ,gL0, VNa, VK, VL, y0);
V = y1(:,1);
m = y1(:,2);
h = y1(:,3);
n = y1(:,4);
tspan2 = [t1(end),T];
y02 = [V(end), m(end),h(end),n(end)];
Iapp2 = 0;
[t2,y2] = ode15s(@HH_completo,tspan2,y02,options, Iapp2, Cm, gNa0, gK0 ,gL0, VNa, VK, VL, y0);
V = [V; y2(:,1)];
m = [m; y2(:,2)];
h = [h; y2(:,3)];
n = [n; y2(:,4)];
t = [t1; t2];

alpha_m = 0.1 * (25 - V) .* (exp((25-V)/(10))-1).^(-1);
alpha_h = 0.07 * exp(-V/20);
alpha_n = 0.01 * (10-V) .* (exp((10-V)/(10))-1).^(-1);
beta_m = 4 * exp(-V/18);
beta_h = (exp((30-V)/(10))+1).^(-1);
beta_n = 0.125 * exp(-V/80);


% nullclines
v = [-10:0.01:100];
N = [0:0.01:1];
for i = 1 : length(v)
        alpha_m = 0.1 * (25 - v(i)) .* (exp((25-v(i))/(10))-1).^(-1);
        alpha_h = 0.07 * exp(-v(i)/20);
        alpha_n = 0.01 * (10-v(i)) .* (exp((10-v(i))/(10))-1).^(-1);
        beta_m = 4 * exp(-v(i)/18);
        beta_h = (exp((30-v(i))/(10))+1).^(-1);
        beta_n = 0.125 * exp(-v(i)/80);
    
        for j = 1:length(N)
        
            dv(i,j)=  -((gK0*N(j).^4.*(v(i)-VK))+(gNa0.*(((alpha_m)./(beta_m+alpha_m)).^3).*(0.8-N(j))*(v(i)-VNa))+(gL0.*(v(i)-VL)));
            dn(i,j)=alpha_n.*(1-N(j))-beta_n.*(N(j));
    end
end


dv_iapp = dv+Iapp;

figure('Name','Nullclines e orbita')
contour(v,N,dv_iapp',[0,0],'b')
hold on
contour(v,N,dv',[0,0],'k')
contour(v,N,dn',[0,0],'r')
plot(V,n)
legend('dv/dt=0', 'dv/dt=0-iapp=0','dn/dt=0','orbita')
title('nullclines e orbita')
xlabel('v')
ylabel('n')

[minV,x_minV] = min(V);

figure
plot(t,V)
hold on
plot(t(x_minV),minV,'*')

%% punto b

clc
clear all 
close all

V0 = -20;
m0 = 5.2934e-2;
h0 = 5.9611e-1;
n0 = 3.1768e-1;

% dati iniziali: 
y0 = [V0,m0,h0,n0];
Iapp = 0; 

Cm = 1;
gNa0 = 120;
gK0 = 36;
gL0 = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;
% 
% T = 150;
% Tstim = 1; % durata dell'impulso di corrente
% tspan = [0,Tstim];
% options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
% [t1,y1] = ode15s(@HH,tspan,y0,options, Iapp, Cm, gNa0, gK0 ,gL0, VNa, VK, VL);
% V = y1(:,1);
% n = y1(:,2);
% tspan2 = [t1(end),T];
% y02 = [V(end), n(end)];
% Iapp2 = 0;
% [t2,y2] = ode15s(@HH,tspan2,y02,options, Iapp2, Cm, gNa0, gK0 ,gL0, VNa, VK, VL);
% V = [V; y2(:,1)];
% n = [n; y2(:,2)];
% t = [t1; t2];



T = 200;
tspan = [0,T];
options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);
[t1,y1] = ode15s(@HH_completo,tspan,y0,options, Iapp, Cm, gNa0, gK0 ,gL0, VNa, VK, VL, y0);
V = y1(:,1);
m = y1(:,2);
h = y1(:,3);
n = y1(:,4);

alpha_m = 0.1 * (25 - V) .* (exp((25-V)/(10))-1).^(-1);
alpha_h = 0.07 * exp(-V/20);
alpha_n = 0.01 * (10-V) .* (exp((10-V)/(10))-1).^(-1);
beta_m = 4 * exp(-V/18);
beta_h = (exp((30-V)/(10))+1).^(-1);
beta_n = 0.125 * exp(-V/80);

% nullclines
v = [-10:0.01:100];
N = [0:0.01:1];
for i = 1 : length(v)
        alpha_m = 0.1 * (25 - v(i)) .* (exp((25-v(i))/(10))-1).^(-1);
        alpha_h = 0.07 * exp(-v(i)/20);
        alpha_n = 0.01 * (10-v(i)) .* (exp((10-v(i))/(10))-1).^(-1);
        beta_m = 4 * exp(-v(i)/18);
        beta_h = (exp((30-v(i))/(10))+1).^(-1);
        beta_n = 0.125 * exp(-v(i)/80);
    
        for j = 1:length(N)
        
            dv(i,j)=  -((gK0*N(j).^4.*(v(i)-VK))+(gNa0.*(((alpha_m)./(beta_m+alpha_m)).^3).*(0.8-N(j))*(v(i)-VNa))+(gL0.*(v(i)-VL)));
            dn(i,j)=alpha_n.*(1-N(j))-beta_n.*(N(j));
    end
end


figure('Name','Nullclines e orbita')
contour(v,N,dv',[0,0],'k')
hold on
contour(v,N,dn',[0,0],'r')
plot(V,n)
legend('dv/dt=0','dn/dt=0','orbita')
title('nullclines e orbita')
xlabel('v')
ylabel('n')

[minV,x_minV] = min(V);

figure
plot(t1,V)
hold on
plot(t1(x_minV),minV,'*')