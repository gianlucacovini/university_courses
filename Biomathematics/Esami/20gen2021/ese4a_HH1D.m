clear
close all
clc

%% ESERCIZIO 4a DEL 20/01/21
%MODELLO HH IN 1D
%PROPAGAZIONE DEL POTENZIALE D'AZIONE LUNGO UNA FIBRA
% Cm*Vt - sigma*Vxx = -(G_Na*m^3*h*(V-V_Na)+G_K*n^4*(V-V_K)+G_L*(V-V_L)) + I_app 
% dm/dt=a_m*(1-m)-b_m*(m)
% dh/dt=a_h*(1-h)-b_h*(h)
% dn/dt=a_n*(1-n)-b_n*(n)

x0=0;
x1=1; %porzione di assone lunga 1cm

%PARAMETRI
sigma=0.01;
Cm=1; % picoFa/cm^2
G_Na=120; %m Ohm^-1/cm^2
G_K=36; %m Ohm^-1/cm^2
G_L=0.3; %m Ohm^-1/cm^2
V_Na=115; %mV
V_K=-12; %mV
V_L=10.6; %mV


%TEMPO
T=20; %ms

%DATI DI NEUMANN (derivate nello spazio al bordo)
BC1=0;
BC2=0;

%MESH SPAZIO-TEMPO
h=0.01;
delta_t=0.05; %Passo temporale da testo
m=T/delta_t;

t=0:delta_t:T; %definisco la mesh in tempo
x=x0:h:x1; %definisco la mesh in spazio

%% Applicare per 1msec uno stimolo Iapp1 = 100ms nei primi 5 nodi dell'assone.
% Descrivere la dinamica generata da questo stimolo e dare una stima della
% velocità di propagazione del fronte generato

%CORRENTE APPLICATA
I_app=zeros(length(x),1);
stimolo=100; %stimolo di durata Iapp1 = 100 mA
Tstim=1; %Applicare per 1ms

%PUNTI INIZIALI
V0=2.7570e-4;
m0=5.2934e-2;
h0=5.9611e-1;
n0=3.1768e-1;

v_k=zeros(length(t),length(x));
v_k(1,:)=V0; 
m_k=zeros(length(t),length(x));
m_k(1,:)=m0;
h_k=zeros(length(t),length(x));
h_k(1,:)=h0;
n_k=zeros(length(t),length(x));
n_k(1,:)=n0;

N=1/h+1; %calcolo il numero di punti
A=(1/h^2)*( 2*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1) );
A(1,2)=-2*(1/h^2);
A(end,end-1)=-2*(1/h^2);
I=eye(N,N);
AA=(Cm*I/delta_t + sigma*A);

%% 
for k=2:length(t)
    
    if t(k)<Tstim
        I_app(find(x<=0.04))=stimolo;       
    else
        I_app(find(x<=0.04))=0;      
    end
    
    
    v_n=v_k(k-1,:)'; % V all'istante precedente
    m_n=m_k(k-1,:)'; % m all'istante precedente
    h_n=h_k(k-1,:)'; % h all'istante precedente
    n_n=n_k(k-1,:)'; % n all'istante precedente
    
    F=(Cm*I/delta_t)*v_n - (G_Na*m_n.^3.*h_n.*(v_n-V_Na)+G_K* n_n.^4.*(v_n-V_K)+G_L*(v_n-V_L)) + I_app;
    F(1)=F(1)+BC1/h^2;
    F(end)=F(end)+BC2/h^2;
    
    v_k(k,:)=AA\F; 
    
    a_m=0.1*(25-v_n).*( exp( (25-v_n)/10 )-1 ).^(-1);
    a_h=0.07*exp(-v_n/20);
    a_n=0.01*(10-v_n).*( exp( (10-v_n)/10 )-1 ).^(-1);

    b_m=4*exp(-v_n/18);
    b_h=(exp((30-v_n)/10)+1).^(-1);
    b_n=0.125*exp(-v_n/80);

    m_k(k,:)= (delta_t*a_m + m_n)./(1 +delta_t*(a_m+b_m)); 
    h_k(k,:)= (delta_t*a_h + h_n)./(1 +delta_t*(a_h+b_h));
    n_k(k,:)= (delta_t*a_n + n_n)./(1 +delta_t*(a_n+b_n));
end  

%% PLOT

figure()
for j=1:length(t)
   plot(x , v_k(j,:) )
   xlim([0,1])
   ylim([-20,140])
   xlabel('x')
   ylabel('V')
   title(['Potenziale nel tempo T=',num2str(t(j))])
   drawnow
end
%Vedo che si genera un potenziale d'azione molto veloce

%Per dare una stima della velocità di propagazione fisso due tempi e vedo a
%che punto dell'assone è il picco massimo del potenziale

%Scelgo tempi così brevi perchè dal grafico del potenziale mi accorgo che a
%7-8 è già concluso

time1=find(t==2); 
time2=find(t==5);

% [massimo1,indice1] = max(v_k(time1,:));
% [massimo2,indice2] = max(v_k(time2,:));
% x1 = x(indice1)
% x2 = x(indice2)
% vel_potenziale = (x2-x1)/(5-2) %[cm/msec]

figure()
subplot(2,1,1)
plot(x , v_k(time1,:) )
xlim([0,1])
ylim([-120,120])
xlabel('x')
ylabel('V')
title(['Potenziale con secondo stimolo nel tempo T=',num2str(t(time1))])

subplot(2,1,2)  
plot(x , v_k(time2,:) )
xlim([0,1])
ylim([-120,120])
xlabel('x')
ylabel('V')
title(['Potenziale con secondo stimolo nel tempo T=',num2str(t(time2))])


[massimo,indice] = max(v_k,[],1);
t_iniziale = t(indice(1));
t_finale = t(indice(end));
vel = 1/(t_finale-t_iniziale)