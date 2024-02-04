function dy = HH_completo(t,y, Iapp, Cm, gNa, gK ,gL, VNa, VK, VL, y0)

alpha_m = 0.1 * (25 - y(1)) .* (exp((25-y(1))/(10))-1)^(-1);
alpha_h = 0.07 * exp(-y(1)/20);
alpha_n = 0.01 * (10-y(1)) .* (exp((10-y(1))/(10))-1)^(-1);
beta_m = 4 * exp(-y(1)/18);
beta_h = (exp((30-y(1))/(10))+1)^(-1);
beta_n = 0.125 * exp(-y(1)/80);


dy(1,1) = 1/Cm * Iapp - 1/Cm * (gNa*y(2)^3*y(3)*(y(1)-VNa)+gK*y(4)^4*(y(1)-VK)+gL*(y(1)-VL));
dy(2,1) = alpha_m * (1-y(2)) - beta_m * y(2); 
dy(3,1) = alpha_h * (1-y(3)) - beta_h * y(3); 
dy(4,1) = alpha_n * (1-y(4)) - beta_n * y(4); 

end 