clear
close all
clc

%% 
x0=0;
x1=1;
%Da testo

T=70;
sigma=0.001;
b=5;
c=1;
beta=0.1;
delta=1;
gamma=0.25;
e=0.1;

%Valori iniziali nel tempo v(x,0) e w(x,0)
v_0=0;
w_0=0;

%Dati di Neumann (derivate nello spazio al bordo)
BC1=0;
BC2=0;

h=1/100; %1/100;
delta_t=0.05; %1/100; 
m=T/delta_t;

t=0:delta_t:T; %definisco la mesh in tempo
x=x0:h:x1; %definisco la mesh in spazio 
%deltax = 0.01 (101 nodi in spazio)

I_app=zeros(length(x),1);
stimolo=2;
Tstim=1;


v_k=zeros(length(t),length(x));
v_k(1,:)=v_0; 
w_k=zeros(length(t),length(x));
w_k(1,:)=w_0;

N=1/h+1; %calcolo il numero di punti
A=(1/h^2)*( 2*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1) );
A(1,2)=-2*(1/h^2);
A(end,end-1)=-2*(1/h^2);
I=eye(N,N);
AA=(I/delta_t + sigma*A);

%%
for k=2:length(t)
    
    if t(k)<=Tstim
        I_app(find(x<=0.04))=stimolo; 
    else
        I_app(find(x<=0.04))=0;
    end
    
    v_n=v_k(k-1,:)'; % V all'istante precedente
    w_n=w_k(k-1,:)'; % W all'istante precedente
    F=(I/delta_t)*v_n + b*v_n.*(v_n-beta).*(delta-v_n) - c*w_n + I_app;
    F(1)=F(1)+BC1/h^2;
    F(end)=F(end)+BC2/h^2;
    
    v_k(k,:)=AA\F;  
    w_k(k,:)= (delta_t/(1+delta_t*e*gamma)) * ( (1/delta_t) * w_n + e * v_n);
    
end               

figure()
for j=1:20:length(t)
   plot(x , v_k(j,:) )
   xlim([0,1])
   ylim([-1,2])
   xlabel('x')
   ylabel('V')
   title(['Potenziale nel tempo T=',num2str(t(j))])
   drawnow
end


[massimo,indice] = max(v_k,[],1);
t_iniziale = t(indice(1));
t_finale = t(indice(end));
vel = 1/(t_finale-t_iniziale);