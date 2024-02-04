%% Es3: HH 0D
clc
close all
clear all

%%  Simulare il PdA:

Cm=1;
gNA=120;
gK=36;
gL=0.3;
Vna=115;
Vk=-12;
VL=10.6;

y0=[2.7570e-04, 5.2934e-02, 5.9611e-01, 3.1768e-01];  % [V0, m0, h0, n0]
Iapp =[0:0.25:250];
tspan=[0:0.05:200];
corrente=0;


vMin=zeros();
vMax=zeros();
for i=1:length(Iapp)

    [t31,y31] = ode15s(@(t,y)odefun(t,y,Iapp(i),Cm, gNA, gK,gL,Vna,Vk,VL),[0:0.05:200],y0);
%     ynew=[y31(end,1),y31(end,2),y31(end,3),y31(end,4)];
%     [t32,y32] = ode15s(@(t,y)odefun3(t,y,corrente,Cm, gNA, gK,gL,Vna,Vk,VL),[0.1:0.05:200],ynew);
%     tot=[t31',t32'];
%     ytot=[y31',y32'];

%     vMin(i)=min(y31(find(t31(150:200)),1));
%     vMax(i)=max(y31(find(t31(150:200)),1));

index = find(t31 == 150);
V_1(i) = min(y31(index:end, 1));
V_2(i) = max(y31(index:end, 1));

figure(1)
plot(t31, y31(:, 1))

end

figure(2)
plot(Iapp,V_1)
grid on
hold on
plot(Iapp,V_2)
xlabel('I_{app}')
ylabel('v(I_{app})')
legend('v_{min}','v_{max}')
sgtitle('Es.3')