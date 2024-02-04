function [dy] = HHkinetik2(t,y,Cm,gNa,gK,gL,VNa,VK,VL,Iapp,Iapp2,Tstim,Tstim2,t_bar)

v = y(1);
n = y(2);

alfa_m = 0.1*(25-v)*(exp((25-v)/10)-1)^-1;
beta_m = 4*exp(-v/18);
minf = alfa_m/(alfa_m+beta_m);
alfa_n = 0.01*(10-v)*(exp((10-v)/10)-1)^-1;
beta_n = 0.125*exp(-v/80);
ninf = alfa_n/(alfa_n+beta_n);
taun = 1/(alfa_n+beta_n);

m = minf;
h = 0.8-n;

INa = gNa*m^3*h*(v-VNa);
IK = gK*n^4*(v-VK);
IL = gL*(v-VL);

if t <= Tstim 
    Iapp = Iapp; 
elseif (t >=t_bar &&  t<=t_bar+Tstim2)
    Iapp = Iapp2;
else
    Iapp = 0;
end


dv = (-(INa+IK+IL)+Iapp)/Cm;
dn = (ninf - n)/taun;

dy = [dv;dn];

end

