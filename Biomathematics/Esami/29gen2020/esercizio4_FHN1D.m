% esercizio 4

%% punto a

clc 
clear all
close all

sigma = 0.001;

alpha = 1.1;

T = 30; % Tmax
m = 600; % numero istanti di tempo

nvect = 99;
 
n = nvect;
u = zeros(n,m+1);
w = zeros(n,m+1);
h = 0.01;
[X]=meshgrid(0:h:1); %creo la griglia ma esce una matrice quindi prendo solo prima riga
mesh=X(1,:);
mesh(end) = [];
mesh(1) = [];

Iapp = zeros(length(mesh),1);
count = 0;
for v=1:length(mesh)
    if mesh(v) >= 0 && mesh(v)<=0.04
       Iapp(v) = alpha;
       count = count + 1;
    else 
        Iapp(v) = 0;
    end
end


a=-1*ones(1,n-1);
c=a;
b=2*ones(1,n);
A=gallery('tridiag',a,b,c);
A(1,2) = -2;
A(end,end-1) = -2;
A = 1/(h^2) * A;

A = sigma * A;
t_index = 1;

beta = 5; 
beta2 = 0.1;
delta2 = 1;
c = 1;
e = 0.1;
gamma = 0.25;

% soluzione al tempo zero al variare di x
u(:,t_index) = zeros(length(mesh),t_index);
w(:,t_index) = zeros(length(mesh),t_index);

I = eye(length(A));


Tstim = 1/4;
% Tstim = 1/2;
% Tstim = 1;
% Tstim = 2;
% Tstim = 4;

t = zeros(1,m+1);
%for negli istanti temporali
for k  = 1 : m
    deltat = T/m;
    
    % calcolo w
     f_w = e *(u(:,k)-gamma*w(:,k));
     U_k_w = w(:,k)+ deltat * f_w;
     w(:,k+1) =  U_k_w;
     
    % calcolo u
    if t(k)<=Tstim
        f= beta * u(:,k) .* (u(:,k)-beta2) .* (delta2-u(:,k)) - c * w(:,k) + Iapp(1); 
        f(count+1:end)= beta * u(count+1:end,k) .* (u(count+1:end,k)-beta2) .* (delta2-u(count+1:end,k)) - c * w(count+1:end,k) + Iapp(end); 
    else
         f= beta * u(:,k) .* (u(:,k)-beta2) .* (delta2-u(:,k)) - c * w(:,k) + Iapp(end);
    end
    
    U_k = u(:,k)/deltat + f;
    A_ =  A + I/deltat;
    u(:,k+1) =  A_\U_k;
    t(k+1) = k * deltat;
    
    
    figure(1)
    plot(mesh,u(:,k+1));
    title("Come varia la V(x) iterando su t");
    xlabel("spazio"), ylabel("V(x,t_I_T_E_R_A_T_O)")
    axis ( [0 1 -0.5 2] )
    pause(0.1)
end


% plot in tempo della soluzione nel punto medio x = 0.5
figure
plot(t,u(51,:))
xlabel('t')
ylabel('u')
title('u vs t per x=0.5')

%% punto b

soglie = [1.1 0.6 0.4 0.2 0.15];
durate = [1/4 1/2 1 2 4];

figure
plot(durate,soglie)
xlabel('durate')
ylabel('soglie')
