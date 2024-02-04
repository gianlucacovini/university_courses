function [u]=Iapp_finale(x,tk)
for j=1:length(x)
    if j >= length(x)-5 && tk<=1
        u(j,1)=-100;
    else
        u(j,1)=0;
end
end
end