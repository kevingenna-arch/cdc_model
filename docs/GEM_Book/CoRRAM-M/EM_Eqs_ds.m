function [f,x] = EM_Eqs_ds();
% Example model from CoRRAM-M user guide. Equations for symbolic differentiation in levels
% Alfred Mauﬂner
% 
% 18 January 2017
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

f=[ y1-(a1^(1-alpha))*(N1^(1-alpha))*(k1^alpha);
    theta*(c1/(1-N1))-(1-alpha)*(y1/N1);    
    a1-astar*exp(z1);
    a1*k2-y1-(1-delta)*k1+c1;
    1-beta*(a1^(-eta))*((c1/c2)^(eta))*(((1-N2)/(1-N1))^(theta*(1-eta)))*(1-delta+alpha*(y2/k2));
    ];
return;
end

