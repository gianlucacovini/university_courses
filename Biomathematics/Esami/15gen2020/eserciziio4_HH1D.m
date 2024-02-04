% esercizio 4

clc 
clear all
% close all

sigma = 0.001;

alpha = 40;

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

m_tempo = 400; % numero istanti di tempo

nvect = 101;
 
u = zeros(nvect,m_tempo+1);
m = zeros(nvect,m_tempo+1);
h = zeros(nvect,m_tempo+1);
n = zeros(nvect,m_tempo+1);
h_mesh = (1/(nvect+1));
[X]=meshgrid(0:h_mesh:1); %creo la griglia ma esce una matrice quindi prendo solo prima riga
mesh=X(1,:);
mesh(end) = [];
mesh(1) = [];

Iapp = [alpha -2*alpha 0];
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


%for negli istanti temporali
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
        f(end-count:end)= 1/Cm * Iapp(2)- 1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(end-count:end,k).^3.*h(end-count:end,k).*(u(end-count:end,k)-VNa)+gK0.*n(end-count:end,k).^4.*(u(end-count:end,k)-VK)+gL0*(u(end-count:end,k)-VL));
    else
         f =  1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
    end
    
    U_k = u(:,k)/deltat + f;
    A_ =  A + I/deltat;
    u(:,k+1) =  A_\U_k;
    t(k+1) = k * deltat;
    
     figure(1)
    plot(mesh,u(:,k+1));
    titoloplot = sprintf('Come varia la V(x) iterando su t\nTempo= %.2f', t(k+1)) ;
    title(titoloplot);
    xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
    axis([0 1 -100 120])
    pause(0.001)
end

Indici=find(tempo~=0);
durata=tempo(Indici(end))-tempo(Indici(1));
velocity=1/durata; %[cm/s];

figure('Name','t=3')
subplot(3,1,1)
plot(mesh,u(:,76))
xlabel('spazio')
ylabel('tempo')
subplot(3,1,2)
plot(mesh,m(:,76))
hold on
plot(mesh,h(:,76))
plot(mesh,n(:,76))
legend('m','h','n')
subplot(3,1,3)
plot(mesh,iNa(:,76))
hold on
plot(mesh,iK(:,76))
plot(mesh,iL(:,76))
plot(mesh,iTot(:,76))
legend('iNa','iK','iL','itot')


figure('Name','t=9')
subplot(3,1,1)
plot(mesh,u(:,226))
xlabel('spazio')
ylabel('tempo')
subplot(3,1,2)
plot(mesh,m(:,226))
hold on
plot(mesh,h(:,226))
plot(mesh,n(:,226))
legend('m','h','n')
subplot(3,1,3)
plot(mesh,iNa(:,226))
hold on
plot(mesh,iK(:,226))
plot(mesh,iL(:,226))
plot(mesh,iTot(:,226))
legend('iNa','iK','iL','itot')


figure('Name','t=11')
subplot(3,1,1)
plot(mesh,u(:,276))
xlabel('spazio')
ylabel('tempo')
subplot(3,1,2)
plot(mesh,m(:,276))
hold on
plot(mesh,h(:,276))
plot(mesh,n(:,276))
legend('m','h','n')
subplot(3,1,3)
plot(mesh,iNa(:,276))
hold on
plot(mesh,iK(:,276))
plot(mesh,iL(:,276))
plot(mesh,iTot(:,276))
legend('iNa','iK','iL','itot')


% Trovare gli istanti di tempo iniziali dei potenziali d'azione

[istanti_provvisori_row , istanti_provvisori_col ] = find(u>=15) ;
% Stimo l'istante di partenza del potenziale a sx
indice_inizio_sx = find(istanti_provvisori_row==1);
indice_tempo_sx = min(istanti_provvisori_col(indice_inizio_sx));
istante_partenza_potenziale_sx = t(indice_tempo_sx) ;
% Stimo l'istante di partenza del potenziale a dx
indice_inizio_dx = find(istanti_provvisori_row==101);
[ indice_tempo_dx , indiceindice_tempo_dx ] = min(istanti_provvisori_col(indice_inizio_dx));
istante_partenza_potenziale_dx = t(indice_tempo_dx) ;

% Se vuole quello del fornte d'onda sviluppato al massimo, allora uso i massimi...
[M,index] = max(u,[],2); % trovo il max e il corrispondente indice per ogni riga della matrice (cioè per ogni x)
istante_partenza_fronte_sx = t(index(1)) ;
istante_partenza_fronte_dx = t(index(end)) ;


% Stimare la velocita di propagazione del fronte;

% Assunzione: supponiamo si muova in maniera uniforme. In realtà andrebbero
% eliminati i transitori iniziale e finale ma volendo fare giusto una stima
% della velocità non è indispensabile e il risultato non cambia di molto.

[M,index] = max(u,[],2); % trovo il max e il corrispondente indice per ogni riga della matrice (cioè per ogni x)


%Considero 50 nodi a dx per stimare la velocità
t_i_sx = t(index(1));
t_f_sx = t(index(50));
delta_t_sx = t_f_sx-t_i_sx;
velocita_fronte_sx = 50*h_mesh/delta_t_sx; %cm/s %Uso 15*h perché è la porzione della fibra che sto considerando per stimare la velocità

%Considero 15 nodi a dx per stimare la velocità
indexx0_8 = find(mesh==0.77) ;

t_i_dx = t(index(end));
t_f_dx = t(index(indexx0_8));
delta_t_dx = abs(t_f_dx-t_i_dx);
velocita_fronte_dx = (length(mesh)-indexx0_8)*h_mesh/delta_t_dx; %cm/s %Uso 15*h perché è la porzione della fibra che sto considerando per stimare la velocità



% indice = find(t==7.9)
% figure('Name','t=8.6')
% subplot(3,1,1)
% plot(mesh,u(:,indice))
% xlabel('spazio')
% ylabel('tempo')
% subplot(3,1,2)
% plot(mesh,m(:,indice))
% hold on
% plot(mesh,h(:,indice))
% plot(mesh,n(:,indice))
% legend('m','h','n')
% subplot(3,1,3)
% plot(mesh,iNa(:,indice))
% hold on
% plot(mesh,iK(:,indice))
% plot(mesh,iL(:,indice))
% plot(mesh,iTot(:,indice))
% legend('iNa','iK','iL','itot')


