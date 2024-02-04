clear
close all

f = @(t, y) -100*y+100*t+1;
tspan = [0 10];
y0 = 3;

y = @(t) y0*exp(1).^(-100*t)+t;

%% AB h=0.0054

h = 0.0054;

[uA1, tA1] = AB3(f, tspan, y0, h);

numItA1 = size(tA1, 2);
for n=1:numItA1
    EA1(n) = uA1(n)-y(tA1(n));
end

%% AB h=0.0055

h = 0.0055;

[uA2, tA2] = AB3(f, tspan, y0, h);

numItA2 = size(tA2, 2);
for n=1:numItA2
    EA2(n) = uA2(n)-y(tA2(n));
end

%% Grafici

figure(1)
subplot(2,1,1)
plot(tA1, uA1, '-o', tA1, y(tA1))
title('Adam-Bashford con h=0.0045')

subplot(2,1,2)
plot(tA2, uA2, '-o', tA2, y(tA2))
title('Adam-Bashford con h=0.0055')