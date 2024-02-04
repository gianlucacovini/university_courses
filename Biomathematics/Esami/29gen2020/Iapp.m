function [u]=Iapp(x,tk, alpha, T_stim)
for j=1:length(x)
    if x(j)<=0.04 && tk<=T_stim % || x(j)>=1-0.04 && tk<=1 % Da decommentare nel caso di doppio stimolo
        u(j,1)=alpha;
    else
        u(j,1)=0;
end
end
end