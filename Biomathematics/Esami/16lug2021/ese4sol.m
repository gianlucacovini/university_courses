%% Esercizio 4 a) - 16 luglio 2021
close all 
clear all

% parametri e dati iniziali
y0=[2.7570e-04, 5.2934e-02, 5.9611e-01, 3.1768e-01]; 

Cm=1; sigma=0.005;
gNa=120; gK=36; gL=0.3; 
VNa=115; VK=-12; VL=10.6;

h=0.01; x=0:h:1;  n=length(x);
T=50;  tau=0.05; t=0:tau:T;   m=length(t);
Iapp = 50;
A=zeros(n,n); % matrice di rigidezza
V=zeros(m,n); M=zeros(m,n); H=zeros(m,n); N=zeros(m,n);
vk=ones(1,length(x))*y0(1); mk=ones(1,length(x))*y0(2); hk=ones(1,length(x))*y0(3); nk=ones(1,length(x))*y0(4);

% matrice A
for i=1:n
    for s=1:n
        if abs(i-s)==0       % per i=j
            A(i,s)=2*sigma*(1/(h^2));  
        elseif abs(i-s)==1   % per i=j+1 oppure i=j-1
            A(i,s)=-1*sigma*(1/(h^2));
        end 
    end 
end 
% A(1,2)=-2*sigma*(1/(h^2));
% A(end,end-1)=-2*sigma*(1/(h^2));
A(1,end)=-1*sigma*(1/(h^2));
A(end,1)=-1*sigma*(1/(h^2));

I=eye(n);
matrice=I+tau*A;

V(1,:)=vk;  M(1,:)=mk;  H(1,:)=hk;  N(1,:)=nk;

for k=1:m-1
   
   f1=(-(gNa*M(k,:).^3.*H(k,:).*(V(k,:)-VNa)+gK*N(k,:).^4.*(V(k,:)-VK)+gL*(V(k,:)-VL)))/Cm; %+Iapp(k,:)
   if t(k)<=5
    f1(18:22)=f1(18:22)+Iapp;
   end 
   
   f2=(0.1*(25-V(k,:)).*(exp((25-V(k,:))/10)-1).^(-1)).*(1-M(k,:))-4*exp(-V(k,:)/18).*M(k,:);
   f3=0.07*exp(-V(k,:)/20).*(1-H(k,:))-(exp((30-V(k,:))/10)+1).^(-1).*H(k,:);
   f4=(0.01*(10-V(k,:)).*(exp((10-V(k,:))/10)-1).^(-1)).*(1-N(k,:))-(0.125*exp(-V(k,:)/80)).*N(k,:);
   
   V(k+1,:)=(matrice)\(V(k,:)'+tau*f1'); % dv/dt
   M(k+1,:)=M(k,:)+tau*f2;               % dm/dt
   H(k+1,:)=H(k,:)+tau*f3;               % dh/dt
   N(k+1,:)=N(k,:)+tau*f4;               % dn/dt
end 

%% Dinamiche
for i= 1:5:length(t)
   figure(1)
   plot(x,V(i,:),'-');
   xlabel("x")
   ylabel("V(x)")
   title("Potenziale d'azione lungo la fibra")
   axis([0 1 -100 120])
end    

gNa_time=gNa*M.^3.*H;
INa=gNa_time.*(V-VNa);
gK_time=gK*N.^4;
IK=gK_time.*(V-VK);
IL=gL*(V-VL);
Iion=INa+IK+IL;

% correnti
for i= 1:5:length(t)
   figure(2)
   subplot 311
   plot(x,INa(i,:),'-');
   xlabel("x")
   ylabel("V(x)")
   xlim([0 1])
   ylim([-1e3 100]) 
   title("I_N_a")
   
   subplot 312
   plot(x,IK(i,:),'-');
   xlabel("x")
   ylabel("V(x)")
   xlim([0 1])
   ylim([0 1e3]) 
   title("I_K")
      
   subplot 313
   plot(x,Iion(i,:),'-');
   xlabel("x")
   ylabel("V(x)")
   xlim([0 1])
   ylim([-300 500])
   title("I_i_o_n")
end   


%% Grafici al tempo t=5.05
figure(3)
% x=0.2
x_def=(find(x==0.2));
subplot 311
plot(t,V(:,x_def),'-');
xlabel("t")
ylabel("V(t)")

subplot 312
plot(t,M(:,x_def),'r-');
hold on
plot(t,H(:,x_def),'g-');
hold on
plot(t,N(:,x_def),'b-');
legend('m','h','n')

subplot 313
plot(t,INa(:,x_def))
hold on
plot(t,IK(:,x_def))
hold on
plot(t,IL(:,x_def))
hold on
plot(t,Iion(:,x_def))
legend('I_N_a','I_K','I_L','I_i_o_n')
sgtitle('x=0.2')

%% Grafici al tempo t=5.05
figure(4)
% t=5.05
t_def=102;
subplot 311
plot(x,V(t_def,:),'-');
xlabel("x")
ylabel("V(x)")
axis([0 1 -20 120])

subplot 312
plot(x,M(t_def,:),'r-');
hold on
plot(x,H(t_def,:),'g-');
hold on
plot(x,N(t_def,:),'b-');
legend('m','h','n')

subplot 313
plot(x,INa(t_def,:))
hold on
plot(x,IK(t_def,:))
hold on
plot(x,IL(t_def,:))
hold on
plot(x,Iion(t_def,:))
legend('I_N_a','I_K','I_L','I_i_o_n')
sgtitle('t=5.05')


%% velocità di propagazione
[massimo,ind]=max(V,[],1); % max e indice corrispondente per ogni colonna: indice dell'istante di tempo in cui c'è stato il picco il quel punto del cavo
t_iniziale=t(ind(1));
t_finale=t(ind(end));
delta_t=t_finale-t_iniziale;
v1=1/delta_t; %cm/msec
