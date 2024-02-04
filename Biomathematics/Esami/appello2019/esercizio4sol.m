% esercizio 4


clc 
clear all
% close all

sigma = 0.001;

alpha = 50;

Cm = 1;

var_gNa0 = [40 60 80 100 120 140];

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

velocity = zeros(1,length(var_gNa0));

for z = 1 : length(var_gNa0)
%for negli istanti temporali
for k  = 1 : m_tempo
    deltat = T/m_tempo;
    
    if max(u(:,k))>70
        tempo(k)=t(k);
    else
        tempo(k)=0;
    end
    
    gNa0 = var_gNa0(z);
    
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
end

Indici=find(tempo~=0);
durata=tempo(Indici(end))-tempo(Indici(1));
velocity(z)=1/durata; %[cm/s];

end

figure
plot(var_gNa0,velocity, '-o')
xlabel('gNa')
ylabel('velocità')


%% punto b


clc 
clear all
close all

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

T = 45;

m_tempo = 900; % numero istanti di tempo
 
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

Iapp = [0 -2*alpha 0];
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

indice_1 = find(mesh>=0.4 & mesh<=0.6);


gNa0_star = 25;


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
        f(indice_1) =  - 1/Cm * (gNa0_star.*m(indice_1,k).^3.*h(indice_1,k).*(u(indice_1,k)-VNa)+gK0.*n(indice_1,k).^4.*(u(indice_1,k)-VK)+gL0.*(u(indice_1,k)-VL));
        f(end-count:end)= 1/Cm * Iapp(2)- 1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(end-count:end,k).^3.*h(end-count:end,k).*(u(end-count:end,k)-VNa)+gK0.*n(end-count:end,k).^4.*(u(end-count:end,k)-VK)+gL0*(u(end-count:end,k)-VL));
    else
         f =  1/Cm * Iapp(end) - 1/Cm * (gNa0.*m(:,k).^3.*h(:,k).*(u(:,k)-VNa)+gK0.*n(:,k).^4.*(u(:,k)-VK)+gL0.*(u(:,k)-VL));
         f(indice_1) =  - 1/Cm * (gNa0_star.*m(indice_1,k).^3.*h(indice_1,k).*(u(indice_1,k)-VNa)+gK0.*n(indice_1,k).^4.*(u(indice_1,k)-VK)+gL0.*(u(indice_1,k)-VL));
    end
    
    U_k = u(:,k)/deltat + f;
    A_ =  A + I/deltat;
    u(:,k+1) =  A_\U_k;
    t(k+1) = k * deltat;
    
     figure(1)
    plot(mesh,u(:,k+1));
    axis([0 1  -20 120])
    title("Come varia la V(x) iterando su t");
    xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
   
%     pause(0.1)
end

Indici=find(tempo~=0);
durata=tempo(Indici(end))-tempo(Indici(1));
velocity=1/durata; %[cm/s];



