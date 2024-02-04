% esercizio4


clc 
clear all
close all

sigma = 0.01;

alpha = 100;

Cm = 1;

gNa0 = 120;

gK0 = 36;
gL0 = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;

V0 = 2.7570e-4;
m0 = 5.2934e-2;
h0 = 5.9611e-1;
n0 = 3.1768e-1;

T = 20;

m_tempo = 500; % numero istanti di tempo

nvect = 99;
 
u = zeros(nvect,m_tempo+1);
m = zeros(nvect,m_tempo+1);
h = zeros(nvect,m_tempo+1);
n = zeros(nvect,m_tempo+1);
h_mesh = (1/(nvect+1));
[X]=meshgrid(0:h_mesh:1); %creo la griglia ma esce una matrice quindi prendo solo prima riga
mesh=X(1,:);
mesh(end) = [];
mesh(1) = [];

Iapp = [alpha 0];
count = 5;


a=-1*ones(1,nvect-1);
c=a;
b=2*ones(1,nvect);
A=gallery('tridiag',a,b,c);
A(1,2) = -2;
A(end,end-1) = -2;
A = 1/(h_mesh^2) * A;

A = sigma * A;
t_index = 1;

% soluzione al tempo zero al variare di x
u(:,t_index) = ones(length(mesh),t_index)*V0;
m(:,t_index) = ones(length(mesh),t_index)*m0;
h(:,t_index) = ones(length(mesh),t_index)*h0;
n(:,t_index) = ones(length(mesh),t_index)*n0;

I = eye(length(A));

t = zeros(1,m_tempo+1);



for k  = 1 : m_tempo
    deltat = T/m_tempo;
    
    if max(u(:,k))>70
        tempo(k)=t(k);
    else
        tempo(k)=0;
    end
    
    
    alpha_m = 0.1 .* (25 - u(:,k)) .* (exp((25-u(:,k))/(10))-1).^(-1);
    alpha_h = 0.07 .* exp(-u(:,k)/20);
    alpha_n = 0.01 .* (10-u(:,k)) .* (exp((10-u(:,k))/(10))-1).^(-1);
    beta_m = 4 .* exp(-u(:,k)/18);
    beta_h = (exp((30-u(:,k))/(10))+1).^(-1);
    beta_n = 0.125 .* exp(-u(:,k)/80);
    
    iNa(:,k)=gNa0*(m(:,k).^3).*h(:,k).*(u(:,k)-VNa);
    iK(:,k)=gK0*(n(:,k).^4).*(u(:,k)-VK);
    iL(:,k)=gL0*(u(:,k)-VL);
    iTot(:,k)=(iNa(:,k)+iK(:,k)+iL(:,k));
    
    % calcolo le variabili di gating 
     f_m = alpha_m .* (1-m(:,k))- beta_m .* m(:,k);
     f_h = alpha_h .* (1-h(:,k))- beta_h .* h(:,k);
     f_n = alpha_n .* (1-n(:,k))- beta_n .* n(:,k);
     
     U_k_m = m(:,k)+ deltat * f_m;
     m(:,k+1) =  U_k_m;
     U_k_h = h(:,k)+ deltat * f_h;
     h(:,k+1) =  U_k_h;
     U_k_n = n(:,k)+ deltat * f_n;
     n(:,k+1) =  U_k_n;
     
    % calcolo u (sarebbe la v)
    if t(k)<=1
        f = 1/Cm * Iapp(1) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
        f(count+1:end)= 1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(count+1:end,k).^3.*h(count+1:end,k).*(u(count+1:end,k)-VNa)+gK0.*n(count+1:end,k).^4.*(u(count+1:end,k)-VK)+gL0*(u(count+1:end,k)-VL));
    else
         f =  1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
    end
    
    U_k = u(:,k)/deltat + f;
    A_ =  A + I/deltat;
    u(:,k+1) =  A_\U_k;
    t(k+1) = k * deltat;
%     
%     figure(1)
%     plot(mesh,u(:,k+1));
%     title("Come varia la V(x) iterando su t");
%     xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
%     pause(0.001)
end

Indici=find(tempo~=0);
durata=tempo(Indici(end))-tempo(Indici(1));
velocity=1/durata; %[cm/s];

[~,indice1] = max(u(1, :));
[~,indice2] = max(u(end, :));
t_in = t(indice1);
t_fin = t(indice2);

vel_potenziale = 1/(t_fin-t_in) %[cm/msec]


% Soluzione nel punto centrale della fibra
indiceX05 = find(mesh==0.5);
VX05 = u(indiceX05,:);
MX05 = m(indiceX05,:);
NX05 = n(indiceX05,:);
HX05 = h(indiceX05,:);

figure("Name","Soluzionen x=0.5");
subplot(3,1,1)
plot(t,VX05);
xlabel("tempo");
ylabel("V(0.5,t)");
title("Soluzionen x=0.5");
subplot(3,1,2)
plot(t,MX05);
hold on
plot(t,NX05);
plot(t,HX05);
xlabel("tempo");
legend('m','n','h')
title("Soluzionen x=0.5");
subplot(3,1,3)
plot(t(1:end-1),iNa(indiceX05,:))
hold on
plot(t(1:end-1),iK(indiceX05,:))
plot(t(1:end-1),iL(indiceX05,:))
plot(t(1:end-1),iTot(indiceX05,:))
legend('iNa','iK','iL','itot')

%% punto b

clc 
clear all
close all

sigma = 0.01;

alpha = 100;

Cm = 1;

gNa0 = 120;

gK0 = 36;
gL0 = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;

V0 = 2.7570e-4;
m0 = 5.2934e-2;
h0 = 5.9611e-1;
n0 = 3.1768e-1;

T = 40;

m_tempo = 800; % numero istanti di tempo

nvect = 99;
 
u = zeros(nvect,m_tempo+1);
m = zeros(nvect,m_tempo+1);
h = zeros(nvect,m_tempo+1);
n = zeros(nvect,m_tempo+1);
h_mesh = (1/(nvect+1));
[X]=meshgrid(0:h_mesh:1); %creo la griglia ma esce una matrice quindi prendo solo prima riga
mesh=X(1,:);
mesh(end) = [];
mesh(1) = [];

Iapp = [alpha 0];
count = 5;


a=-1*ones(1,nvect-1);
c=a;
b=2*ones(1,nvect);
A=gallery('tridiag',a,b,c);
A(1,2) = -2;
A(end,end-1) = -2;
A = 1/(h_mesh^2) * A;

A = sigma * A;
t_index = 1;

% soluzione al tempo zero al variare di x
u(:,t_index) = ones(length(mesh),t_index)*V0;
m(:,t_index) = ones(length(mesh),t_index)*m0;
h(:,t_index) = ones(length(mesh),t_index)*h0;
n(:,t_index) = ones(length(mesh),t_index)*n0;

I = eye(length(A));

t = zeros(1,m_tempo+1);

Iapp2 = 450;

for k  = 1 : m_tempo
    deltat = T/m_tempo;
    
    if max(u(:,k))>70
        tempo(k)=t(k);
    else
        tempo(k)=0;
    end
    
    
    alpha_m = 0.1 .* (25 - u(:,k)) .* (exp((25-u(:,k))/(10))-1).^(-1);
    alpha_h = 0.07 .* exp(-u(:,k)/20);
    alpha_n = 0.01 .* (10-u(:,k)) .* (exp((10-u(:,k))/(10))-1).^(-1);
    beta_m = 4 .* exp(-u(:,k)/18);
    beta_h = (exp((30-u(:,k))/(10))+1).^(-1);
    beta_n = 0.125 .* exp(-u(:,k)/80);
    
    iNa(:,k)=gNa0*(m(:,k).^3).*h(:,k).*(u(:,k)-VNa);
    iK(:,k)=gK0*(n(:,k).^4).*(u(:,k)-VK);
    iL(:,k)=gL0*(u(:,k)-VL);
    iTot(:,k)=(iNa(:,k)+iK(:,k)+iL(:,k));
    
    % calcolo le variabili di gating 
     f_m = alpha_m .* (1-m(:,k))- beta_m .* m(:,k);
     f_h = alpha_h .* (1-h(:,k))- beta_h .* h(:,k);
     f_n = alpha_n .* (1-n(:,k))- beta_n .* n(:,k);
     
     U_k_m = m(:,k)+ deltat * f_m;
     m(:,k+1) =  U_k_m;
     U_k_h = h(:,k)+ deltat * f_h;
     h(:,k+1) =  U_k_h;
     U_k_n = n(:,k)+ deltat * f_n;
     n(:,k+1) =  U_k_n;
     
    % calcolo u (sarebbe la v)
    if t(k)<=1
        f = 1/Cm * Iapp(1) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
        f(count+1:end)= 1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(count+1:end,k).^3.*h(count+1:end,k).*(u(count+1:end,k)-VNa)+gK0.*n(count+1:end,k).^4.*(u(count+1:end,k)-VK)+gL0*(u(count+1:end,k)-VL));
    else
         f =  1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
    end
    
    if 8 <= t(k) && t(k) <= 10
        f(48:52)= 1/Cm * Iapp2 - 1/Cm * (gNa0.*m(48:52,k).^3.*h(48:52,k).*(u(48:52,k)-VNa)+gK0.*n(48:52,k).^4.*(u(48:52,k)-VK)+gL0*(u(48:52,k)-VL));
    end
    
    U_k = u(:,k)/deltat + f;
    A_ =  A + I/deltat;
    u(:,k+1) =  A_\U_k;
    t(k+1) = k * deltat;
    
%     figure(1)
%     plot(mesh,u(:,k+1));
%     title("Come varia la V(x) iterando su t");
%     xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
    
end
for k  = 1 : 2 : m_tempo
   figure(1)
    plot(mesh,u(:,k+1));
    title("Come varia la V(x) iterando su t");
    xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
    axis([0 1 -30 120])
    pause(0.001)
    end
Indici=find(tempo~=0);
durata=tempo(Indici(end))-tempo(Indici(1));
velocity=1/durata; %[cm/s];




% Soluzione nel punto centrale della fibra
indiceX05 = find(mesh==0.5);
VX05 = u(indiceX05,:);
MX05 = m(indiceX05,:);
NX05 = n(indiceX05,:);
HX05 = h(indiceX05,:);


figure("Name","Soluzionen x=0.5");
subplot(3,1,1)
plot(t,VX05);
xlabel("tempo");
ylabel("V(0.5,t)");
title("Soluzionen x=0.5");
subplot(3,1,2)
plot(t,MX05);
hold on
plot(t,NX05);
plot(t,HX05);
xlabel("tempo");
legend('m','n','h')
title("Soluzionen x=0.5");
subplot(3,1,3)
plot(t(1:end-1),iNa(indiceX05,:))
hold on
plot(t(1:end-1),iK(indiceX05,:))
plot(t(1:end-1),iL(indiceX05,:))
plot(t(1:end-1),iTot(indiceX05,:))
legend('iNa','iK','iL','itot')



indice_t = find(t==7.5);
figure("Name","Soluzionen t=7.5");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')

indice_t = find(t==8);
figure("Name","Soluzionen t=8");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')


indice_t = find(t==9);
figure("Name","Soluzionen t=9");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')


indice_t = find(t==9.25);
figure("Name","Soluzionen t=9.2");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')


indice_t = find(t==9.5);
figure("Name","Soluzione in t=9.5");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')

indice_t = find(t==9.65);
figure("Name","Soluzionen t=9.65");
subplot(2,1,1)
plot(mesh,u(:,indice_t));
xlabel("spazio");
ylabel("V");
subplot(2,1,2)
plot(mesh,m(:,indice_t));
hold on
plot(mesh,n(:,indice_t));
plot(mesh,h(:,indice_t));
xlabel("spazio");
legend('m','n','h')