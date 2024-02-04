clear
close all
clc

%% variabili e condizioni iniziali - 1
s0 = 5;
i0 = 2;
e0 = 1;
k1 = 102;
k_1 = 50;
k2 = 1;
k3 = 26;
k_3 = 50;
ki=k_3/k3;
ks=k_1/k1;

e = @(t,y) e0-y(3)-y(4)-y(5); %per non doverlo sostituire a mano

% in ordine s i x y z
f = @(t,y) [-k1.*e(t,y).*y(1)+k_1.*y(3)+k_1.*y(5)-k1.*y(4).*y(1);...
            -k3.*e(t,y).*y(2)+k_3.*y(4)-k3.*y(3).*y(2)+k_3.*y(5);...
            -k_1.*y(3)+k1.*e(t,y).*y(1)-k3.*y(3).*y(2)+k_3.*y(5)-k2.*y(3);...
            k3.*e(t,y).*y(2)-k_3.*y(4)+k_1.*y(5)-k1.*y(4).*y(1);...
            -k_3.*y(5)+k3.*y(3).*y(2)-k_1.*y(5)+k1.*y(4).*y(1)];

options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);

[T,S]=ode15s(f,[0 100], [s0;i0;0;0;0], options);

%% plot soluzioni - 2

figure(1)
plot(T,S)
title("Inibizione allosterica")
legend("s","i","x","y","z",Location="best")

%% plot velocità - 3

V = ((k2*e0)./(1+S(:,2)./ki)).*(S(:,1)./(S(:,1)+ks));
V=V';

figure(2)
plot(S(:,1),V)
axis([0 5 0 inf])
grid on
title("Velocità in funzione di s")

figure(3)
plot3(S(:,1),S(:,2),V)
xlabel("s")
ylabel("i")
zlabel("V")
grid on
title("Velocità in funzione di s e i")





