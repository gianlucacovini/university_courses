function [V, m, h, n, X,T]=HH1D_solver(f,I,IT,n_x,n_t,v_x0,v_x1,v0,m0,h0,n0,Sig,C_m,alpha_m,beta_m,alpha_h,beta_h,alpha_n,beta_n)
% I è l'intervallo di spazio
% T è l'intervallo di tempo
% t,h sono il passo in tempo e in spazio
% u_x0=u_x(0,t), u_x1=u_x(1,t), u0=u(x,0)

t=(IT(2)-IT(1))/n_t;
delta_x=(I(2)-I(1))/(n_x+1); % ex h

X=linspace(I(1),I(2),n_x); % n_x nodi totali
T=linspace(IT(1),IT(2),n_t); % n_t nodi totali

X=X'; T=T';
X=[X(1)-delta_x;X;X(end)+delta_x]; % aggiungiamo i ghost points

V_old=v0(X); %primo passo nel tempo
m_old=m0(X);
h_old=h0(X);
n_old=n0(X);
%fare un tensore se li si vuole salvare al passare del tempo

for k=1:n_t % è il ciclo sul tempo

    F=f(V_old, m_old, h_old, n_old); % F viene calcolata ad ogni tempo
    F(1)=F(1)-2*v_x0*delta_x; % Controlla il segno -
    F(n_x+2)=F(n_x+2)+2*v_x1*delta_x;
    F=F+Iapp(X,T(k));

    A=diag(2*ones(1,n_x+2)) + diag(-1*ones(1,n_x+1),1) + diag(-1*ones(1,n_x+1),-1);
    A(1,2)=-2;
    A(n_x+2,n_x+1)=-2;
    A=(Sig./delta_x.^2).*A;
  
    %il sistema: (I+tA)U_k+1=U_k+tF_k
    M=eye(n_x+2,n_x+2);
    V_new=(C_m*M+t.*A)\(C_m*V_old+t.*F); % +b*V_old.*(V_old-beta).*(delta-V_old)-c*W_old));

        V(:,k)=V_old;

    V_old=V_new; % aggiorno il dato

    m_new = (M+t.*alpha_m(V_old)+t.*beta_m(V_old))\(M*m_old+t.*alpha_m(V_old));

        m(:, k) = m_old;

    m_old = m_new;

    h_new = (M+t.*alpha_h(V_old)+t.*beta_h(V_old))\(M*h_old+t.*alpha_h(V_old));

        h(:, k) = h_old;

    h_old = h_new;

    n_new = (M+t.*alpha_n(V_old)+t.*beta_n(V_old))\(M*n_old+t.*alpha_n(V_old));

        n(:, k) = n_old;

    n_old = n_new;

end



end