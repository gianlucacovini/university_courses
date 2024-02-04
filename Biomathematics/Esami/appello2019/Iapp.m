function [u]=Iapp(x,tk)
for j=1:length(x)
    if x(j)<=0.04 && tk<=1
        u(j,1)=50;
    else
        u(j,1)=0;
end
end
end