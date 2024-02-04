function [u]=Iapp(x,tk)
for j=1:length(x)
    if j<=5 && tk<=1
        u(j,1)=100;
    elseif j>=48 && j<=52 && tk>=7 && tk<=9
        u(j, 1) = 450;
    else
        u(j,1)=0;
end
end
end