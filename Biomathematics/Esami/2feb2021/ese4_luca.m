clc
clear all
close all

sigma = 0.001;
bi = 5;
ci = 1;
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

I = 2;

n = 79;
tau = 0.02;
T = 30;
nt = T/tau;
v0 = 0;
w0 = 0;

h = 1/(n+1);
x = [0:h:0.8];
t = [0:tau:T];
tindex = 1;
tstim=1;

Iapp = zeros(length(x),length(t));
v = zeros(length(x),length(t));

%costruisco matrice Iapp
for i=1:length(x)
    for j=1:length(t)
        if i<=4 && t(j)<=tstim

            Iapp(i,j) = I;
        end
        v(i,j)=x(i)/t(j);
    end
end

%costruisco matrice A
A = zeros(n+2,n+2); %n+2 = ndof = end

for i=1:(n+2)
    for j=1:(n+2)
        if i==j          
          A(i,j) = -2;
        elseif abs(i-j)==1
          A(i,j) = +1;
        end
    end   
end
A(1,2)=+2;
A(n+2,n+1)=+2;
% A(1,n+2)=1;
% A(n+2,1)=1;
v_a = zeros(length(x),length(t));
v_a(:,1) = ones(length(x),1)*v0;
w_a = zeros(length(x),length(t));
w_a(:,1) = ones(length(x),1)*w0;

B = eye(n+2) - tau*sigma*A/h^2;%come esercitazione 5
b=zeros(length(x),1); 
c=zeros(length(x),1);

for kk=1:nt  
    vj=v_a(:,tindex);
    wj=w_a(:,tindex);
    t_noto_v = (bi*vj.*(vj-beta).*(delta-vj)-ci*wj); % è la f(uj) dell'equazione
    t_noto_w = e*(vj-gamma*wj);
    b = (vj + tau*t_noto_v + tau*Iapp(:,tindex+1));
    c = (wj + tau*t_noto_w);
    tindex=tindex+1;
    v_a(:,tindex) = B\b;
    w_a(:,tindex) = c;
end

% figure()
% a=t==15;
% plot(x,v_a(:,a))
% title('v nello spazio a t=15ms')
% figure()
% plot(t,v_a(51,:))
% title('v nel tempo al nodo medio 51')


for i=1:20:length(t)
    figure(3)
  
    plot(x',v_a(:,i));
    xlabel('X')
    ylabel('Vappr')
    ylim([-0.5 3.5])
    xlim([0 1])
    title('Grafico della Vappr')

end
% con Iapp=1 non parte il potenziale e potrebbe essere per il codice che è
% un eulero esplicito invece che implicito 