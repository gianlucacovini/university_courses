clear
close all

% Dovrebbe essere fatto (al limite verifica la funzione)

%% FHN1D

I=[0 1];
IT=[0 20];

n=101; 
m=400; % Per avere passo delta_t = 0.05. Forse va preso 401...

sig = 0.001; 
b = 5; 
c = 1; 
beta = 0.1; 
delta = 1; 
gamma = 0.25; 
e = 0.1;


v0=@(x) zeros(size(x));
w0=@(x) zeros(size(x));

f=@(v, w) b.*v.*(v-beta).*(delta-v)-c*w; % dato di FHN1D
v_x0=0; v_x1=0; %dato al bordo di Neumann

T_stim = 1/4;

for alpha=0.9:0.005:0.92
    [V,W,X,T]=FHN1D_solver(f,I,IT,n,m,v_x0,v_x1,v0,w0,sig, gamma, e, alpha, T_stim);

    figure(1)
    hold on
    txt = ['I_{app} =',num2str(alpha)];
    plot(T,V(floor(n/2), :),"DisplayName",txt)
    ylim([-0.3 1.3])
    title("potenziale al variare di I_{app}")
    xlabel("t")
    ylabel("V")
end
hold off
legend show

%% Risultati

% Per 1/4 è 0.915
% Per 1/2 è 0.475
% Per 1 è 0.255
% Per 2 è 0.145
% Per 4 è 0.095

%% Plot risultati

Ts = [1/4 1/2 1 2 4];
Is = [1 0.5 0.26 0.15 0.095];

figure(2)
loglog(Ts, Is, Ts, 0.25*1./Ts)

%% Risultati

% La curva è molto simile a 1/(4*T_stim)
