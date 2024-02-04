function [g]=G_Na(x, gNaSalto)
for j=1:length(x)
    if x(j) > 0.4 && x(j) <= 0.6
        g(j)=gNaSalto;
    else
        g(j) = 120;
    end
end
end