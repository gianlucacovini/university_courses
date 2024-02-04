function [T,V,m,h,n]= HH_solver(F,I,v0,m0,h0,n0,T_di)

options = odeset('RelTol',5.e-13,'AbsTol',[1.e-13 1.e-13 1.e-13 1.e-13],"InitialStep",1.e-5,"MaxStep",5);

l=floor(I(2)/T_di); %quanti problemi dobbiamo risolvere per fare passi di lunghezza T_di
%l'ultimo passo potrebbe sforare, andiamo per difetto
a=I(1);b=I(2);
clear I
I=zeros(1,l+1);
I(1)=a;
for i=2:l+1
I(i)=I(i-1)+T_di;
end
I(l+1)=b; %nel caso ci fosse dell'eccesso l'ultimo intervallino più lungo

k=1;

v_k=v0;
m_k=m0;
h_k=h0;
n_k=n0;

V=[];
m=[];
h=[];
n=[];
T=[];



%% servono per aggiornare F, per ora li metto qui


%%
%osserviamo che alla fine F non dipende dal tempo se non in Iapp, allora
%possiamno semplicemente risolvere il sistema su [0 I(k+1)-I(k)] che è
%intervallo di lunghezza uguale ma ha il vantaggio di tenere Iapp classica
%e di non dover variare F ad ogni iterazione

while k<=l

[T_k,X_k] = ode15s(F,[0 I(k+1)-I(k)],[v_k,m_k,h_k,n_k],options);

%bisogna solo sistemare in modo che l'intervallo [I(k) I(k+1)] alla fine
%abbia i nodi disposti come quelli di T_k, basta fare:

T_k=T_k+I(k);


%attacco il pezzo di soluzione
V=[V;X_k(:,1)];
m=[m;X_k(:,2)];
h=[h;X_k(:,3)];
n=[n;X_k(:,4)];
T=[T;T_k];

v_k=V(end);
m_k=m(end);
h_k=h(end);
n_k=n(end);

k=k+1;
end
