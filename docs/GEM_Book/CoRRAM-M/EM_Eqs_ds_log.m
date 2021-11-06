function [f,x] = EM_Eqs_ds_log();
% Example model from CoRRAM-M user guide. Difference stationary growth. Equations for symbolic differentiation in logs
% Alfred Mauﬂner
% 
% 02 Feburary 2017
%

% Parameters
alpha=sym('alpha');
beta=sym('beta');
eta=sym('eta');
delta=sym('delta');
theta=sym('theta');
astar=sym('astar');

% variables of the model, 1 refers to period t and 2 to period t+1 variables
syms k1 z1 a1 y1 c1 N1 k2 z2 a2 y2 c2 N2;

x=[k2 z2 a2 y2 c2 N2 k1 z1 a1 y1 c1 N1];

% equations of the model

f=[ exp(y1)-exp((1-alpha)*a1+(1-alpha)*N1+alpha*k1);
    theta*exp(c1)-(1-alpha)*exp(y1-N1)*(1-exp(N1));    
    exp(a1)-astar*exp(z1);
    exp(a1+k2)-exp(y1)-(1-delta)*exp(k1)+exp(c1);
    1-beta*exp(-eta*a1)*exp(eta*(c1-c2))*(((1-exp(N2))/(1-exp(N1)))^(theta*(1-eta)))*(1-delta+alpha*exp(y2-k2));
    ];
return;
end

