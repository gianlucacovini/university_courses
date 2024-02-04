close all
clear
clc

%% variabili e condizioni iniziali - 1
s0 = 5;
i0 = 1;
e0 = 1;
k1 = 102;
k_1 = 50;
k2 = 1;
k3 = 26;
k_3 = 50;

%y(1)=s, y(2)=i, y(3)=c1, y(4)=c2
f = @(t,y) [-k1.*y(1).*(e0-y(3)-y(4))+k_1.*y(3) ;...
    -k3.*y(2).*(e0-y(3)-y(4))+k_3.*y(4) ;...
    k1.*y(1).*(e0-y(3)-y(4))-(k_1+k2).*y(3) ;...
    k3.*y(2).*(e0-y(3)-y(4))-k_3.*y(4)];

options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);

[T,S]=ode15s(f,[0 100], [s0;i0;0;0], options);

%% spprossimazione quasi stazionaria - 2

km=(k_1+k2)/k1;
ki=k_3/k3;

%y(1)=2; y(2)=i
c1 = @(y) ((e0.*y(:,1))./(y(:,1)+km.*(1+y(:,2)./ki)));
c2 = @(y) ((e0.*y(:,2))./(y(:,2)+ki.*(1+y(:,1)./km)));

%sostituendo a mano c1 e c2 per evitare problemi
g = @(t,y) [-k1.*y(1).*(e0-((e0.*y(1))./(y(1)+km.*(1+y(2)./ki)))-((e0.*y(2))./(y(2)+ki.*(1+y(1)./km))))+k_1.*((e0.*y(1))./(y(1)+km.*(1+y(2)./ki))) ;...
    -k3.*y(2).*(e0-((e0.*y(1))./(y(1)+km.*(1+y(2)./ki)))-((e0.*y(2))./(y(2)+ki.*(1+y(1)./km))))+k_3.*((e0.*y(2))./(y(2)+ki.*(1+y(1)./km)))] ;
options1 = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);
[T,Ysi]=ode15s(g,T, [s0;i0], options1);

Y(:,1)=Ysi(:,1);
Y(:,2)=Ysi(:,2);
Y(:,3)=c1(Ysi(:,1:2));
Y(:,4)=c2(Ysi(:,1:2));


%% plot - 3
figure(1)
plot(T,S,T,Y,"--")
title("soluzione esatta e quasi-stazionaria")
legend("s","i","c_{1}","c_{2}","s_{u}","i_{u}","c_{1}_{u}","c_{2}_{u}")

%% plot dell'errore - 4
E=abs(S-Y);

figure(2)
plot(T,E);
title("Grafico errore assoluto")
legend("err_{s}","err_{i}","err_{c_{1}}","err_{c_{2}}")
%% grafico della velocità V - 5

V(:) = k2.*S(:,3); %k2*c1
V=V';
Vu(:) = ((k2*e0)./(1+Y(:,2)./(k1))).*(Y(:,1)./(Y(:,1)+km.*(1+Y(:,2)./ki)));
Vu=Vu';

figure(3)
subplot(2,1,1)
plot3(S(:,1),S(:,2),V)
grid on
hold on
plot3(Y(:,1),Y(:,2),Vu)
hold off
title("Grafico delle velocità")
legend("V","V_{u}",Location="best")

%% grafico errore V - 6
Ev(:)=abs(V(:)-Vu(:));

subplot(2,1,2)
plot(Y(:,1),Ev)
axis([0 5 -.1 .1])
title("Grafico errori della velocità")
legend("V-V_{u}",Location="best")

%% senza inibitore - 7
i0=0;
%risolvo con il nuovo dato
[T,S_i]=ode15s(f,T, [s0;i0;0;0], options);

figure(4)

subplot(1,3,1)
plot(T,S)
title("Con inibitore")
legend("s","i","c_{1}","c_{2}")

subplot(1,3,2)
plot(T,S_i)
title("Senza inibitore")
legend("s_{i}","i_{i}","c_{1_{i}}","c_{2_{i}}")


V_i(:) = k2.*S_i(:,3); %k2*c1
V_i=V_i';

subplot(1,3,3)
plot(Y(:,1),V,Y(:,1),V_i)
axis([0 5 0 inf])
title("Confronto velocità")
legend("V","V_{i}")

sgtitle("Confronto con e senza inibitore")