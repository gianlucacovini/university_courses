close all
%% dati
C_m = 1;
G_na = 120;
G_k = 36;
G_l = 0.3;
V_na = 115;
V_k = -12;
V_l = 10.6;

a_m=@(v) 0.1.*(25-v)./(exp((25-v)./10)-1);
a_h=@(v) 0.07.*exp(-v./20);
a_n=@(v) 0.01.*(10-v)./(exp((10-v)./10)-1);
b_m=@(v) 4.*exp(-v./18);
b_h=@(v) 1./(exp((30-v)./10)+1);
b_n=@(v) 0.125.*exp(-v./80);

I=[0 111]; % Per risolvere il problema dell'approssimazione per difetto
% che tagliava le ultime ripetizioni allunghiamo il tempo su cui calcola
% del necessario e tagliamo poi l'asse x nel grafico

%% quarta parte, vario il numero di applicazioni di Iapp
%Per farlo chiamo una function esterna che risolve k problemi consecutivi in cui il
%primo millisecondo è attiva Iapp e nel resto non lo è


T_di=[20;17;15;12;10;8;7;6;5;4;3;2;1.5];

figure(1)
for i=1:length(T_di)
    
    Iapp=@(t) 100*(t<0.1);
    
    F=@(t,x) [-G_na.*x(3).*(x(2).^3).*(x(1)-V_na)-G_k.*(x(4).^4).*(x(1)-V_k)-G_l.*(x(1)-V_l)+Iapp(t);...
                a_m(x(1)).*(1-x(2))-b_m(x(1))*x(2);...
                a_h(x(1)).*(1-x(3))-b_h(x(1))*x(3);...
                a_n(x(1)).*(1-x(4))-b_n(x(1))*x(4)];
    
    v0=2.7570e-4; 
    m0=5.2934e-2;
    h0=5.9611e-1;
    n0=3.1768e-1;
    
    
    [T,V,m,h,n]= HH_solver(F,I,v0,m0,h0,n0,T_di(i));
    
    %% altri valori per i plot
    
    g_na=G_na.*(m.^3).*h;
    g_k=G_k.*n.^4;
    
    I_k=g_k.*(V-V_k);
    I_na=g_na.*(V-V_na);
    I_l=G_l.*(V-V_l);
    I_ion=I_na+I_k+I_l;
    
    Tau_n= 1./(a_n(V)+b_n(V));
    Tau_m= 1./(a_m(V)+b_m(V));
    Tau_h= 1./(a_h(V)+b_h(V));
    
    n_inf= a_n(V)./(a_n(V)+b_n(V));
    m_inf= a_m(V)./(a_m(V)+b_m(V));
    h_inf= a_h(V)./(a_h(V)+b_h(V));
    
    
    %% plot
    subplot(4, 4, i)
    txt = ['#T_{DI} =',num2str(T_di(i))];
    plot(T,V)
    xlim([0 100])
    title(txt)
    xlabel("t")
    ylabel("V")
end

%% OSSERVAZIONI

% Effetto refrattarietà: se si dà un secondo impulso quando c'è del
% potenziale d'azione non se ne genera un altro, a volte non si genera
% nulla. Il sistema è già eccitato non si può fare di più
