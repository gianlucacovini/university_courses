clear
close all
clc

%% dati
C_m = 1;
G_na = 120;
G_k = 36;
G_l = 0.3;
V_na = 115;
V_k = -12;
V_l = 10.6;

count = 1;

Iapp_int = 0:0.25:250;


a_m=@(v) 0.1.*(25-v)./(exp((25-v)./10)-1);
a_h=@(v) 0.07.*exp(-v./20);
a_n=@(v) 0.01.*(10-v)./(exp((10-v)./10)-1);
b_m=@(v) 4.*exp(-v./18);
b_h=@(v) 1./(exp((30-v)./10)+1);
b_n=@(v) 0.125.*exp(-v./80);

v0=2.7570e-04;
m0=5.2934e-2;
h0=5.9611e-1;
n0=3.1768e-1;

for Iapp=Iapp_int

%% variabili V,m,h,n --> fai function

F=@(t,x) [-G_na.*x(3).*(x(2).^3).*(x(1)-V_na)-G_k.*(x(4).^4).*(x(1)-V_k)-G_l.*(x(1)-V_l)+Iapp;...
            a_m(x(1)).*(1-x(2))-b_m(x(1))*x(2);...
            a_h(x(1)).*(1-x(3))-b_h(x(1))*x(3);...
            a_n(x(1)).*(1-x(4))-b_n(x(1))*x(4)];


% options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[T,X] = ode15s(F,0:0.05:200,[v0;m0;h0;n0]);

%% altri valori per i plot
V=X(:,1);

index = find(T == 150);
V_1(count) = min(V(index:end));
V_2(count) = max(V(index:end));

count = count + 1;

end
hold off

%% plot
figure(1)
plot(Iapp_int, V_1, Iapp_int, V_2)