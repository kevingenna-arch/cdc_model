function fx = EM_Sys_ds_log(w,eq_n,Par)
% Example model from CoRRAM-M user guide. Equations for numeric differentiation, difference stationary version, log
% variables

%{
    Copyright Alfred Maußner
 
    30 January 2017

     Input:     w,  2*(nx+ny+nz) by 1 vector, the stationary solution of the model
              eq_n, scalar, if >0, return equation no eq_n.
               Par, 5 by 1 vector, values of the model's parameter (see below)

     Output:  fx, nx+ny by 1 vector, the left-hand side of the model's equations evaluated at the stationary solution

%}

%Parameters
astar=Par(1);
alpha=Par(2);
beta=Par(3);
delta=Par(4);
eta=Par(5);
theta=Par(6);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
w=exp(w);

k2=w(1);    k1=w(7);
z2=w(2);    z1=w(8);
a2=w(3);    a1=w(9);
y2=w(4);    y1=w(10);
c2=w(5);    c1=w(11);
N2=w(6);    N1=w(12);


% equations of the model
fx=ones(5,1);

fx(1)=y1-(a1^(1-alpha))*(N1^(1-alpha))*(k1^alpha);
fx(2)=theta*(c1/(1-N1))-(1-alpha)*(y1/N1);
fx(3)=a1-astar*z1;
fx(4)=a1*k2-y1-(1-delta)*k1+c1;
fx(5)=1-beta*(a1^(-eta))*((c1/c2)^(eta))*(((1-N2)/(1-N1))^(theta*(1-eta)))*(1-delta+alpha*(y2/k2));

if eq_n>0; fx=fx(eq_n); end;

return;
end

