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

Iapp = 0;

a_m=@(v) 0.1.*(25-v)./(exp((25-v)./10)-1);
a_h=@(v) 0.07.*exp(-v./20);
a_n=@(v) 0.01.*(10-v)./(exp((10-v)./10)-1);
b_m=@(v) 4.*exp(-v./18);
b_h=@(v) 1./(exp((30-v)./10)+1);
b_n=@(v) 0.125.*exp(-v./80);

I=[0 20];

%% seconda parte, variare il valore di v0 per ottenere potenziale d'azione

F=@(t,x) [-G_na.*x(3).*(x(2).^3).*(x(1)-V_na)-G_k.*(x(4).^4).*(x(1)-V_k)-G_l.*(x(1)-V_l)+Iapp;...
            a_m(x(1)).*(1-x(2))-b_m(x(1))*x(2);...
            a_h(x(1)).*(1-x(3))-b_h(x(1))*x(3);...
            a_n(x(1)).*(1-x(4))-b_n(x(1))*x(4)];

v0=4:1:8; %trovare valore che attiva PA
m0=5.2934e-2;
h0=5.9611e-1;
n0=3.1768e-1;

for i=1:length(v0)
    
    options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
    [T,X] = ode15s(F,I,[v0(i);m0;h0;n0],options);
    
    V=X(:,1);
    m=X(:,2);
    h=X(:,3);
    n=X(:,4);
    
    %% altri valori per i plot
    
    g_na=G_na.*(m.^3).*h;
    g_k=G_k.*n.^4;
    
    I_k=g_k.*(V-V_k);
    I_na=g_na.*(V-V_na);
    I_l=G_l.*(V-V_l);
    I_ion=I_na+I_k+I_l;
    
    Tau_n= 1./(a_n(V)+b_n(V));
    Tau_m= 1./(a_m(V)+b_m(V));
    Tau_h= 1./(a_h(V)+b_h(V));
    
    n_inf= a_n(V)./(a_n(V)+b_n(V));
    m_inf= a_m(V)./(a_m(V)+b_m(V));
    h_inf= a_h(V)./(a_h(V)+b_h(V));
    
    
    %% plot
    figure(1)
    hold on
    txt = ['V0 =',num2str(v0(i))];
    plot(T,V,"DisplayName",txt)
    title("potenziale al variare di V_{0}")
    xlabel("t")
    ylabel("V")
end
hold off
legend show
