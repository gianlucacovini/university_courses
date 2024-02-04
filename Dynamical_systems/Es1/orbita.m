function []=orbita(a,t)
    tspan = 0:0.01:t;
    u1 = @(t) exp(1).^(-t)+sin(t);
    u2 = @(t) exp(1).^(-t) + cos(a.*t);
    
    figure
    plot(u1(tspan),u2(tspan))
    title('Orbita u_1 contro u_2')
    xlabel('u_1(t)')
    ylabel('u_2(t)')
end

% La funzione non è periodica ma quasi periodica perché al passare del
% tempo la funzione si sposta, non torna sempre sugli stessi punti ma si
% sposta all'esterno. Dalle sole componenti sembrava periodica la funzione.
% Questo indipendentemente dal parametro.