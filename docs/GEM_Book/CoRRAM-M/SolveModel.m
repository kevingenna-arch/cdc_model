function [Solution,rc]=SolveModel(Model,Par)
% Provides first through third order perturbation solution of DSGE model
%{
    Copyright Alfred Mauﬂner

    First version: 11 Januar 2017 (developed from SolveModelA3)

    Purpose:

    Solve solve for the linear, quadratic, and third-order coefficients of a 
    perturbation solution of a DSGE model. Depending on the flag Model.numeric
    the derivatives of the system of equations are either approximated numerically
    (works for Model.order=1 or =2 but not for Model.order=3) or computed from
    computed from the analytic Jacobian, Hessian, and third-order derivative matrices
    of the model's equations.

    The program recognizes the partition of the vector of states into a vector of
    enodgenous states x(t) and purely exogenous shocks z(t).

    The model is given as

           Eg(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))=0
           z(t+1)-Rho*z(t)-Omega*epsilon(t+1)=0

    The program expects that the ordering of the variables in the matrices of first-, second-, and third-order
    derivatives of g is  (x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))=0.

    Depending on the setting of Model.numeric the vector valued function g must be defined analytically 
    in the file Model.Equations or in the function handle Model.equations.
 
    The structure Par stores the symbols used for the model's parameters and their related numeric values.
    The structure Model.Var stores - among other information used to plot impulse responses and
    to print simulation results - the symbolic names of the model's variables and their respective
    numeric values in the model's stationary solution.


    Let w(t)=[x(t);z(t)]. The approixmate solution for the vector x(t+1) is given by

           x(t+1)-x=Hx_w*[w(t)-w]                                           (linear part)
                    + 0.5*Hx_ss + kron(eye(nx),[w(t)-w])*Hx_ww*[w(t)-w]     (quadratic part)
                    +(1/6)*Hx_sss + kron(eye(nx),[w(t)-w]*Hx_ssx            (cubic part)
                    + 0.5*kron(kron(eye(nx),[w(t)-w]),[w(t)-w])*Hx_www*[w(t)-w]

    The solution for the vector y(t) has the same structure:

           y(t)  -y=Hy_w*[w(t)-w]                                           (linear part)
                    + 0.5*Hy_ss + kron(eye(ny),[w(t)-w])*Hy_ww*[w(t)-w]     (quadratic part)
                    +(1/6)*Hy_sss + kron(eye(ny),[w(t)-w]*Hy_ssx            (cubic part)
                    + 0.5*kron(kron(eye(ny),[w(t)-w]),[w(t)-w])*Hy_www*[w(t)-w]

     Input:    Model, an instance of the DSGE class
               Model.nx: integer, the number of endogenous state variables (elements in x(t))
               Model.ny: integer, the number of control or jump variables (elements in y(t))
               Model.nz: integer, the number of shocks (elements in z(t))
               Model.nu: integer, the number of static equations
              Model.Rho: nz by nz matrix with the autoregressive coefficients of the VAR(1) defining the dynamics of the shocks
            Model.Omega: nz by nz matrix so that Sigma=E(Omega'*epsilon')(Omega*epsilon)  is the covariance matrix of shocks Omega*epsilon
             Model.Skew: nz by nz^2 matrix (skewness of disturbances)
              Model.Var: structure of length nx+ny+nz which stores information about the model's variables,
               the members of the structure needed here are
               . Symbol
               . Star
            
          Par: structure whose length dependes on the number of parameters that appear in
               defining the model's equations, both of its members are needed here:
                . Symbol
                . Value
             
        Model.Equations: string, the name (except for the .m extension) of the file which holds the model's equations
                         or file handel so that Model.Equations([x;z;y;x;z;y]) returns Eg(.) at
                         the stationary solution.

     Output: Solution, a structur with fields:

               Hx_w:     nx by nw matrix
               Hy_w:     ny by nw matrix
               Hx_ww:    nx*nw by nw matrix
               Hy_ww:    ny*nw by nw matrix
               Hx_ss:    nx by 1 matrix
               Hy_ss:    ny by 1 matrix
               Hx_www:   nx*nw^2 by nw matrix
               Hy_www:   ny*nw^2 by nw matrix
               Hx_sss:   nx by 1 matrix
               Hy_sss:   ny by 1 matrix
               Hx_ssw:   nx*nw by 1 matrix
               Hy_ssw:   ny*nw by 1 matrix

               rc:       integer, if equal to 0 a solution has been found, else the index of the error message in Messages
               Messages: cell array with error messages, see code right below


%}

Solution=struct('Hx_w',0,'Hy_w',0,'Hx_ww',0,'Hy_ww',0,'Hx_ss',0,'Hy_ss',0,...
    'Hx_www',0,'Hy_www',0,'Hx_ssw',0,'Hy_ssw',0,'Hx_sss',0,'Hy_sss',0);

% get handle of LogFile.txt
%Model.lfh=fopen('LogFile.txt','a+');
Model.lfh=fopen(strcat(Model.outfile,'_LogFile.txt'),'w');
lfh=Model.lfh;
fprintf(lfh,strcat(char(datetime('now')),'\n'));

% initialize the error code
rc=0;

% check if Par is present
if nargin==2; present=true; else present=false; end;

%{
 Compute the Jacobian matrix, the Hessian matrix, and if required, also the matrix of third-order derivatives.
 After the execution of the following script the numeric Jacobian D is available. If Model.order=2, also the numeric Hessian H is available.
 If Model.order=3, also the numeric matrix of third-order derivatives T is available.
%}
nx=Model.nx;
ny=Model.ny;
nz=Model.nz;
nu=Model.nu;
nw=nx+nz;
nxy=nx+ny;
na=nxy+nz;
Rho=Model.Rho;

GetDHT;

if rc~=0; fclose(lfh); return; end;

% Check for invalid Jacobian
if all(all(isnan(D))) ||   all(all(isinf(D)));
    fprintf(lfh,'Invalid Jacobian\n');
    flclose(lfh);
    rc=3;
    return;
end;

% Obtain linear part of the solution
if Model.reduced;
    [Hx_w,Hy_w,rc]=Linear1(Model,D);
else
    [Hx_w,Hy_w,rc]=Linear2(Model,D);
end;
if rc~=0;
    fclose(lfh);
    return;
end;

Solution.Hx_w=Hx_w;
Solution.Hy_w=Hy_w;

if Model.order==1;    
    fclose(lfh);
    return;
end;

% Obtain quadratric part
Omega=Model.Omega;
Muss=Model.Muss;
Syl=Model.syl;

Quadratic;

Solution.Hx_ww=Hx_ww;
Solution.Hy_ww=Hy_ww;
Solution.Hx_ss=Hx_ss;
Solution.Hy_ss=Hy_ss;

if Model.order<3;   
    fclose(lfh);
    return;
end;

% Obtain cubic part
Skew=Model.Skew;
Cubic;    
Solution.Hx_www=Hx_www;
Solution.Hy_www=Hy_www;
Solution.Hx_ssw=Hx_ssw;
Solution.Hy_ssw=Hy_ssw;
Solution.Hx_sss=Hx_sss;
Solution.Hy_sss=Hy_sss;
fclose(lfh);

return;

end
