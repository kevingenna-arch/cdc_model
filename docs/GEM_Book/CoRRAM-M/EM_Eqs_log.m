function [f,x] = EM_Eqs_log();
% Example model from CoRRAM-M user guide. Equations for symbolic differentiation in logs
% Alfred Mauﬂner
% 
% 13 January 2017
%

% Parameters
alpha=sym('alpha');
beta=sym('beta');
eta=sym('eta');
delta=sym('delta');
theta=sym('theta');

% variables of the model, 1 refers to period t and 2 to period t+1 variables
syms k1 z1 y1 c1 n1 k2 z2 y2 c2 n2;

x=[k2 z2 y2 c2 n2 k1 z1 y1 c1 n1];

% equations of the model

f=[ exp(y1)-exp(z1+n1*(1-alpha)+alpha*k1);
    theta*(exp(c1)/(1-exp(n1)))-(1-alpha)*exp(y1-n1);
    exp(k2)-exp(y1)-(1-delta)*exp(k1)+exp(c1);
    1-beta*exp(eta*(c1-c2))*(((1-exp(N2))/(1-exp(N1)))^(theta*(1-eta)))*(1-delta+alpha*exp(y2-k2));
    ];
return;
end

