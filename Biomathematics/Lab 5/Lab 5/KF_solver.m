function [U,X,T]=KF_solver(f,I,IT,n,m,u_x0,u_x1,u0,Sig,alpha)
% I è l'intervallo di spazio
% IT è l'intervallo di tempo
% t,h sono il passo in tempo e in spazio
% u_x0=u_x(0,t), u_x1=u_x(1,t), u0=u(x,0)

t=(IT(2)-IT(1))/m;
h=(I(2)-I(1))/(n+1); % da controllare se n+1

X=linspace(I(1),I(2),n); % n nodi totali
T=linspace(IT(1),IT(2),m); % m nodi totali

X=X'; T=T';
X=[X(1)-h;X;X(end)+h]; % aggiungiamo i ghost points

U_old=u0(X); %primo passo nel tempo
%fare un tensore se li si vuole salvare al passare del tempo

for k=1:m % è il ciclo sul tempo

    F=f(U_old); % F viene calcolata ad ogni tempo
    F(1)=F(1)-2*u_x0*h; % Correzione di F(1) per condizioni di Neumann
    F(n+2)=F(n+2)+2*u_x1*h; % Correzione dell'ultimo punto (ghost point) per condizioni di Neumann
    F=F+Iapp(X,T(k),alpha);

    A=diag(2*ones(1,n+2)) + diag(-1*ones(1,n+1),1) + diag(-1*ones(1,n+1),-1);
    A(1,2)=-2; % Correzione ghost points
    A(n+2,n+1)=-2;
    A=(Sig./h.^2).*A;
  
    %il sistema: (I+tA)U_{k+1}=U_k+tF_k
    M=eye(n+2,n+2);
    U_new=(M+t.*A)\(U_old+t.*F);

        U(:,k)=U_old; % Salvo le soluzioni prima di sovrascriverle

    U_old=U_new; % aggiorno il dato

end



end