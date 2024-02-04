% Non sono per niente sicuro, non riesco a capire...

function [V, W, X,T]=FHN1D_solver_periodic(f,I,T,h,t,v_x0,v_x1,v0,w0,Sig,b, c, beta, delta, gamma, e)
% I è l'intervallo di spazio
% T è l'intervallo di tempo
% t,h sono il passo in tempo e in spazio
% u_x0=u_x(0,t), u_x1=u_x(1,t), u0=u(x,0)

X=I(1):h:I(2); % n+1 nodi totali
T=T(1):t:T(2); % m nodi totali

X=X'; T=T';

V_old = v0(X); %primo passo nel tempo
W_old = w0(X);
%fare un tensore se li si vuole salvare al passare del tempo
m = length(T);
n = length(X);

for k=1:m % è il ciclo sul tempo

    F=f(V_old, W_old); % F viene calcolata ad ogni tempo
    F=F+Iapp(X,T(k));

    A=diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);
    A(1,end)=-1;
    A(end,1)=-1;
    A=(Sig./h.^2).*A;
  
    %il sistema: (I+tA)U_k+1=U_k+tF_k
    M=eye(n, n);
    V_new=(M+t.*A)\(V_old+t.*F); % +b*V_old.*(V_old-beta).*(delta-V_old)-c*W_old));

        V(:,k)=V_old;

    V_old=V_new; % aggiorno il dato

    W_new = (M+t.*e*gamma.*M)\(W_old+e*t.*V_old);

        W(:,k) = W_old;

    W_old = W_new;

end



end