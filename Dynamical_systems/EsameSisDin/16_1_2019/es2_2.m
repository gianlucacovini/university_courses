clear
close all

tspan = [0 10];

eps = 1;

for k=1:10
    y0 = [(2-(k-1))/5 0 0];
    
    odefun = @(t, y) [(-5*y(1)-y(1)*y(2)+5*y(2)^2+y(3))/(eps+y(2)*y(3)-y(1));...
        (10*y(1)-y(1)*y(2)-10*y(2)^2+y(3))/(eps-y(2)*y(3)+y(1));...
        (y(1)*y(2)-y(3))/(eps-y(2)*y(3)+y(1))];
    
    [t, y] = ode15s(odefun, tspan, y0);
    plot3(y(:,1), y(:,2), y(:,3))
    grid on
    title('Traiettorie nello spazio delle fasi')
    xlabel('y(1)')
    ylabel('y(2)')
    zlabel('y(3)')
    hold on
end