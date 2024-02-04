function [V, W, X,T]=FHN1D_solver(f,I,T,n,m,v_x0,v_x1,v0,w0,Sig, gamma, e, alpha, T_stim)
% I è l'intervallo di spazio
% T è l'intervallo di tempo
% t,h sono il passo in tempo e in spazio
% u_x0=u_x(0,t), u_x1=u_x(1,t), u0=u(x,0)

t=(T(2)-T(1))/m;
h=(I(2)-I(1))/(n+1); % da controllare se n+1

X=linspace(I(1),I(2),n); % n nodi totali
T=linspace(T(1),T(2),m); % m nodi totali

X=X'; T=T';
X=[X(1)-h;X;X(end)+h]; % aggiungiamo i ghost points

V_old=v0(X); %primo passo nel tempo
W_old = w0(X);
%fare un tensore se li si vuole salvare al passare del tempo

for k=1:m % è il ciclo sul tempo

    F=f(V_old, W_old); % F viene calcolata ad ogni tempo
    F(1)=F(1)-2*v_x0*h;
    F(n+2)=F(n+2)+2*v_x1*h;
    F=F+Iapp(X,T(k),alpha, T_stim);

    A=diag(2*ones(1,n+2)) + diag(-1*ones(1,n+1),1) + diag(-1*ones(1,n+1),-1);
    A(1,2)=-2;
    A(n+2,n+1)=-2;
    A=(Sig./h.^2).*A;
  
    %il sistema: (I+tA)U_k+1=U_k+tF_k
    M=eye(n+2,n+2);
    V_new=(M+t.*A)\(V_old+t.*F); % +b*V_old.*(V_old-beta).*(delta-V_old)-c*W_old));

        V(:,k)=V_old;

    V_old=V_new; % aggiorno il dato

    W_new = (M+t.*e*gamma.*M)\(W_old+e*t.*V_old);

        W(:,k) = W_old;

    W_old = W_new;

end



end