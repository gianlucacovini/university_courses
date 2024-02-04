clear
close all

tspan = [0 20];

odefun = @(t, y) [-y(1)*(1-y(2)-y(3))+10^y(3)*(y(2)+y(3))-10^5*y(1)*y(2);...
    y(1)*(1-y(2)-y(3))-10^y(3)*y(2)-10^y(3)*y(1)*y(2)+102*y(3);...
    10^5*y(1)*y(2)-102*y(3)];

for k=1:5
    y0 = [2*k 0 0]';
    
    [t, y] = ode15s(odefun, tspan, y0);

    N(k) = size(t,1);
    
    plot3(y(:,1), y(:,2), y(:,3))
    hold on
    title('Traiettorie')
    xlabel('t')
    ylabel('y')
end