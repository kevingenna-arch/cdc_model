function fx = EM_Sys(w,eq_n,Par)
% Example model from CoRRAM-M user guide. Equations for numeric differentiation

%{
    Copyright Alfred Maußner
 
    13 January 2017

     Input:     w,  2*(nx+ny+nz) by 1 vector, the stationary solution of the model
              eq_n, scalar, if >0, return equation no eq_n.
               Par, 5 by 1 vector, values of the model's parameter (see below)

     Output:  fx, nx+ny by 1 vector, the left-hand side of the model's equations evaluated at the stationary solution

%}

alpha=Par(1);
beta=Par(2);
delta=Par(3);
eta=Par(4);
theta=Par(5);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
K2=w(1);    K1=w(6);
z2=w(2);    z1=w(7);
Y2=w(3);    Y1=w(8);
C2=w(4);    C1=w(9);
N2=w(5);    N1=w(10);

% equations of the model
fx=ones(4,1);

fx(1)=Y1-exp(z1)*(N1^(1-alpha))*(K1^alpha);
fx(2)=theta*(C1/(1-N1))-(1-alpha)*(Y1/N1);
fx(3)=K2-Y1-(1-delta)*K1+C1;
fx(4)=1-beta*(C1/C2)^(eta)*(((1-N2)/(1-N1))^(theta*(1-eta)))*(1-delta+alpha*(Y2/K2));

if eq_n>0; fx=fx(eq_n); end;

return;
end

