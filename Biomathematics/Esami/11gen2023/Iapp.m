function [u]=Iapp(x,tk)
for j=1:length(x)
    if j>=6 && j<=10 && tk<=4
        u(j,1)=50;
    else
        u(j,1)=0;
end
end
end