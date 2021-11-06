function [f,x] = EM_Eqs();
% Example model from CoRRAM-M user guide. Equations for symbolic differentiation in levels
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
syms K1 z1 Y1 C1 N1 K2 z2 Y2 C2 N2;

x=[K2 z2 Y2 C2 N2 K1 z1 Y1 C1 N1];

% equations of the model

f=[ Y1-exp(z1)*(N1^(1-alpha))*(K1^alpha);
    theta*(C1/(1-N1))-(1-alpha)*(Y1/N1);
    K2-Y1-(1-delta)*K1+C1;
    1-beta*((C1/C2)^(eta))*(((1-N2)/(1-N1))^(theta*(1-eta)))*(1-delta+alpha*(Y2/K2));
    ];
return;
end

