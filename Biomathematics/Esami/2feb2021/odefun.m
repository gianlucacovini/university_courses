function dydt = odefun(t,y,Iapp,Cm, gNA, gK,gL,Vna,Vk,VL)
    
    dydt = zeros(4,1);
    dydt(1) =-(gNA*(y(2)^3)*y(3)*(y(1)-Vna)+gK*y(4)^4*(y(1)-Vk)+gL*(y(1)-VL))/Cm + Iapp/Cm;
    
    alfam= 0.1*(25-y(1))*(exp((25-y(1))/10)-1)^-1;
    betam= 4*exp(-y(1)/18);
    dydt(2) = alfam*(1-y(2)) - betam*y(2);
    
    alfah=0.07*exp(-y(1)/20);
    betah=(exp((30-y(1))/10) + 1)^-1;
    dydt(3)=alfah*(1-y(3))-betah*y(3);
    
    alfan=0.01*(10-y(1))*(exp((10-y(1))/10)-1)^-1;
    betan=0.125*exp(-y(1)/80);
    dydt(4)=alfan*(1-y(4))-betan*y(4);

end