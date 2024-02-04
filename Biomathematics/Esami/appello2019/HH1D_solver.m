function [V_vec, m_vec, h_vec, n_vec, X, T]=HH1D_solver(f,I,IT,h,delta_t,v_x0,v_x1,V0,m0,h0,n0,sigma,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n)

x0=I(1);
x1=I(2); 

%MESH SPAZIO-TEMPO
T_fin = IT(2);
m=T_fin/delta_t;

T=0:delta_t:T_fin; %definisco la mesh in tempo
X=x0:h:x1; %definisco la mesh in spazio

V_vec=zeros(length(T),length(X));
V_vec(1,:)=V0(X); 
m_vec=zeros(length(T),length(X));
m_vec(1,:)=m0(X);
h_vec=zeros(length(T),length(X));
h_vec(1,:)=h0(X);
n_vec=zeros(length(T),length(X));
n_vec(1,:)=n0(X);

n=1/h+1; %calcolo il numero di punti
A=(1/h^2)*(2*diag(ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1));
A(1,2)=-2*(1/h^2);
A(end,end-1)=-2*(1/h^2);
I=eye(n,n);
AA=(C_m*I/delta_t + sigma*A);

%% 
for k=2:length(T)
    
    V_old=V_vec(k-1,:)'; % V all'istante precedente
    m_old=m_vec(k-1,:)'; % m all'istante precedente
    h_old=h_vec(k-1,:)'; % h all'istante precedente
    n_old=n_vec(k-1,:)'; % n all'istante precedente
    
    F=(C_m*I/delta_t)*V_old + f(V_old, m_old, h_old, n_old) + Iapp(X, T(k-1));
    F(1)=F(1)+v_x0/h^2;
    F(end)=F(end)+v_x1/h^2;
    
    V_vec(k,:)=AA\F; 

    m_vec(k,:)= (delta_t*alpha_m(V_old) + m_old)./(1 +delta_t*(alpha_m(V_old)+beta_m(V_old))); 
    h_vec(k,:)= (delta_t*alpha_h(V_old) + h_old)./(1 +delta_t*(alpha_h(V_old)+beta_h(V_old)));
    n_vec(k,:)= (delta_t*alpha_n(V_old) + n_old)./(1 +delta_t*(alpha_n(V_old)+beta_n(V_old)));
end  