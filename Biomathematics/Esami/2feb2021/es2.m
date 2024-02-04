clc
clear all
close all

v = linspace(-80,80,17);
Ik = [0.0 0.1 0.4 1.9 7.1 24.7 74.1 172.5 294.2 394.3 466.4 522.9 572.8 620.2 666.7 712.9 759.0];

Eresting = -85;

gk = Ik./(v-Eresting);

subplot(1,2,1)
plot(v,Ik)
axis([-80 80 0 800])
subplot(1,2,2)
plot(v,gk)
axis([-80 80 0 5])

Gk = max(gk);
I = Ik(find(v==-10));

n = (I/(Gk*(-10+85))^1/3);