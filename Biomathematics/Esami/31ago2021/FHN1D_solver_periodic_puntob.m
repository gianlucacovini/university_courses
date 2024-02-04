function [V_vec, w_vec, X, T]=FHN1D_solver2_periodic_puntob(f,I,IT,h,delta_t,V0,w0,sigma,gamma,e)

x0=I(1);
x1=I(2); 

%MESH SPAZIO-TEMPO
T_fin = IT(2);
m=T_fin/delta_t;

T=0:delta_t:T_fin; %definisco la mesh in tempo
X=x0:h:x1; %definisco la mesh in spazio

V_vec=zeros(length(T),length(X));
V_vec(1,:)=V0(X); 
w_vec=zeros(length(T),length(X));
w_vec(1,:)=w0(X);
w_vec(1,13:18)= 0.5; 

n=length(X); 
A=(1/h^2)*(2*diag(ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1));
A(1,end)=-1*(1/h^2);
A(end,1)=-1*(1/h^2);
I=eye(n,n);
AA=(I/delta_t + sigma*A);

%% 
for k=2:length(T)
    
    V_old=V_vec(k-1,:)'; % V all'istante precedente
    w_old=w_vec(k-1,:)'; % w all'istante precedente
    
    F=(I/delta_t)*V_old + f(V_old, w_old) + Iapp(X, T(k-1)); 
    
    V_vec(k,:)=AA\F; 
    
    if k <= 5
        w_vec(k,:)= (I+delta_t.*e*gamma.*I)\(w_old+e*delta_t.*V_old); 
        w_vec(k,13:18)= 0.5; 
    end
end  

end