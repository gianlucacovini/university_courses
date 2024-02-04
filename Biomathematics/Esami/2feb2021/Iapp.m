function [u]=Iapp(x,tk)
for j=1:length(x)
    if j>=18 && j<=22 && tk<=1
        u(j,1)=2;
    else
        u(j,1)=0;
end
end
end