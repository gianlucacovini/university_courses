%% esercizio 3

clear all
% close all
clc


% Lo 0 del potenziale equivale ai -70mV della cellula, a riposo, mentre il
% massimo è 40mV nella cellula. La soglia per far scattare un potenziale è
% a -55mV.
% Quindi, nel caso del modello, il potenziale va da 0 a circa 110, ovvero
% dal valore di riposo (-70) al massimo (+40mV). Il gap è di 110mV, quindi
% andrò da 0 a 110 circa, e la treshold per considerare un potenziale è di
% -55, ovvero di 15mV nel modello.


I_app_soglia = 10 ; %mA
Tstim = 1 ; %ms

% Per la I applicata non parte nessun Potenziale d'Azione



sigma = 0.01;
beta = 4; 
beta2 = 0.1;
delta2 = 1;
c = 1;
e = 0.1;
gamma = 0.25;


n = 99;
tau = 0.05;
T = 100;
nt = T/tau;
h = 1/(n+1);
x = 0:h:1;
t = 0:tau:T;


Iapp = zeros(length(x),length(t));

for j=1:length(t)
    if t(j)<= Tstim
        Iapp ( 18 : 22 , j ) = I_app_soglia;
    end
end


% Costruisco la matrice A, mettendo le condizioni al bordo periodiche

A = zeros(n+2,n+2);
for i=1:(n+2)
    for j=1:(n+2)
        if i==j          
          A(i,j) = -2;
        elseif abs(i-j)==1
          A(i,j) = +1;
        end
    end   
end
% A(1,2)=+2;
% A(n+2,n+1)=+2;

A(end,1)= 1 ;
A(1,end)= 1 ;

% La matrice A finale sarebbe la mia B
B = eye(n+2) /tau - sigma*A/h^2;

% Condizioni iniziali

u = zeros( length(x) , length(t) );
w = zeros( length(x) , length(t) );




for kk = 1 : nt  
    
     % calcolo w com eulero semi-ESPLICITO
     w( : , kk+1 ) = w(:,kk) + tau * e *( u( : , kk ) - gamma * w(:,kk));
    
    % calcolo u
    b = u(:,kk) / tau + beta * u(:,kk) .* (u(:,kk)-beta2) .* (delta2-u(:,kk)) - c * w(:,kk) + Iapp( : , kk+1);
    u( : , kk+1 ) = B\b;

end

periodoplot = 1 ;
% 
for i=1: periodoplot :length(t)    
    figure(1)
    plot(x',u(:,i),'-ro');
    titolotempo = sprintf('Propagazione di 1 PdA\n Tempo: %.2f', t(i)) ;
    title(titolotempo)
    xlabel("Assone, spazio (cm)"), ylabel("V(x,t")
    axis ( [0 1 -0.5 1.5] )
    pause(0.01)
end

% Plot in spazio e tempo del potenziale
it15 = find(t==15) ;
figure(1)
plot(t,u(51,:));
titolotempo = sprintf('PdA nel punto medio al variare del Tempo') ;
title(titolotempo)
xlabel("Tempo ms"), ylabel("V(t)")
% axis ( [0 1 -0.5 1.5] )


figure(2)
plot(x',u(:,t(it15)));
titolotempo = sprintf('PdA lungo l"assone al Tempo: %.2f', t(it15)) ;
title(titolotempo)
xlabel("Assone, spazio (cm)"), ylabel("V(x)")
% axis ( [0 1 -0.5 1.5] )

% La corrente applicata non è sufficiente a generare un Potenziale d'azione
% lungo la fibra del neurone. Infatti il valore massimo del potenziale
% risulta molto più piccolo rispetto a quello normalizzato che si ottiene
% con il modello FHN con un potenziale d'azione attivo.


%% Punto b)
% b) Ripetere il punto a) forzando la variabile di recovery al valore w = 0.5 nei nodi 13 - 18 per
% la durata dello stimolo. Spiegare la dinamica ottenuta

clear all
% close all
clc


% Lo 0 del potenziale equivale ai -70mV della cellula, a riposo, mentre il
% massimo è 40mV nella cellula. La soglia per far scattare un potenziale è
% a -55mV.
% Quindi, nel caso del modello, il potenziale va da 0 a circa 110, ovvero
% dal valore di riposo (-70) al massimo (+40mV). Il gap è di 110mV, quindi
% andrò da 0 a 110 circa, e la treshold per considerare un potenziale è di
% -55, ovvero di 15mV nel modello.


I_app_soglia = 1 ; %mA
Tstim = 1 ; %ms

% Per la I applicata non parte nessun Potenziale d'Azione



sigma = 0.01;
beta = 4; 
beta2 = 0.1;
delta2 = 1;
c = 1;
e = 0.1;
gamma = 0.25;


n = 99;
tau = 0.05;
T = 100;
nt = T/tau;
h = 1/(n+1);
x = 0:h:1;
t = 0:tau:T;


Iapp = zeros(length(x),length(t));


% Condizioni iniziali

u = zeros( length(x) , length(t) );
w = zeros( length(x) , length(t) );


for j=1:length(t)
    if t(j)<= Tstim
        Iapp ( 18 : 22 , j ) = I_app_soglia;
    end
end


% Costruisco la matrice A, mettendo le condizioni al bordo periodiche

A = zeros(n+2,n+2);
for i=1:(n+2)
    for j=1:(n+2)
        if i==j          
          A(i,j) = -2;
        elseif abs(i-j)==1
          A(i,j) = +1;
        end
    end   
end
% A(1,2)=+2;
% A(n+2,n+1)=+2;

A(end,1)= 1 ;
A(1,end)= 1 ;

% La matrice A finale sarebbe la mia B
B = eye(n+2) /tau - sigma*A/h^2;




for kk = 1 : nt      
     
    if t(kk)<= Tstim
        w ( 13 : 18 , kk ) = 0.5 ;
    end
    % calcolo w com eulero semi-ESPLICITO
    w( : , kk+1 ) = w(:,kk) + tau * e *( u( : , kk ) - gamma * w(:,kk));
    

     
    % calcolo u
    b = u(:,kk) / tau + beta * u(:,kk) .* (u(:,kk)-beta2) .* (delta2-u(:,kk)) - c * w(:,kk) + Iapp( : , kk+1);
    u( : , kk+1 ) = B\b;
end
% 
% periodoplot = 1 ;
% 
% for i=1: periodoplot :length(t)    
%     figure(7)
%     plot(x',u(:,i),'-ro');
%     titolotempo = sprintf('Propagazione di 1 PdA\n Tempo: %.2f', t(i)) ;
%     title(titolotempo)
%     xlabel("Assone, spazio (cm)"), ylabel("V(x,t")
%     axis ( [0 1 -0.5 1.5] )
%     pause(0.01)
% end


% Plot in spazio e tempo del potenziale
it15 = find(t==15) ;
figure(3)
plot(t,u(51,:));
titolotempo = sprintf('PdA nel punto medio al variare del Tempo') ;
title(titolotempo)
xlabel("Tempo ms"), ylabel("V(x)")
% axis ( [0 1 -0.5 1.5] )


figure(4)
plot(x',u(:,t(it15)));
titolotempo = sprintf('PdA lungo l"assone al Tempo: %.2f', t(it15)) ;
title(titolotempo)
xlabel("Assone, spazio (cm)"), ylabel("V(t)")
% axis ( [0 1 -0.5 1.5] )


% Il fatto che la variabile w di resting sia posta a un valore così alto per
% 1 ms nei nodi precedenti rispetto a quelli in cui si applica la corrente
% impedisce ulteriormente la generazione del potenziale d'azione sulla
% fibra. Infatti, come si può notare dai grafici il valore massimo
% raggiunto dal potenziale è molto inferiore rispetto al punto a.
% Inoltre se si analizza la matrice u si nota che i valori del poteniale
% negli istanti di tempo di applicazione della corrente nei nodi dove w è
% già attiva, sono di ben 2 ordini di grandezza inferiori.


