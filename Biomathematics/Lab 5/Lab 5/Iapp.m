function [u]=Iapp(x,tk,alpha)
for j=1:length(x)
    if x(j)>=0 && x(j)<=0.04 && tk<=1
        u(j,1)=alpha;
    else
        u(j,1)=0;
end
end
end