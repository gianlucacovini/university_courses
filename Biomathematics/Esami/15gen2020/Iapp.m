function [u]=Iapp(x,tk)
for j=1:length(x)
    if j<=5 && tk<=1
        u(j,1)=40;
    elseif j>=length(x)-5 && tk<=1
        u(j,1) = -80;
    else
        u(j,1)=0;
    end
end
end