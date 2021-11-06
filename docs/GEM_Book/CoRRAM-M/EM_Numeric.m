% EM_Numeric.m
%
%{
 
    Author:         Alfred Maußner
    First Version:  17 January 2017
    Revision:       30 Januar 2017 (option for difference stationary growth and log specification added) 

    Purpose:

    This Matlab script tests and solves the example model from the CoRRAM-M users guide using
    numerical differentiation.

%}

% options:
ds=false;  % if true, difference stationary growth model
log=false; % if true, the model will be formulated in logs

% Calibration
a=0.004;
alpha=0.3;
beta=0.99;
delta=0.025;
eta=2.0;
theta=2.07;
if ds;
    rho=0;
    sigma=0.018;
else
    rho=0.95;
    sigma=0.007;
end;

% Stationary solution
if ds;
    astar=exp(a);
    yk=(astar^eta-beta*(1-delta))/(alpha*beta);
    ck=yk-(astar-1+delta);
    temp=(theta/(1-alpha))*(ck/yk);
    Nstar=1/(1+temp);
    kstar=astar*(yk^(1/(alpha-1)))*Nstar;
    ystar=yk*kstar;
    cstar=ck*kstar;
    %theta=(yk/ck)*((1-Nstar)/Nstar)*(1-alpha);
else
    YK=(1-beta*(1-delta))/(alpha*beta);
    CK=YK-delta;
    temp=(theta/(1-alpha))*(CK/YK);
    Nstar=1/(1+temp);
    Kstar=(YK^(1/(alpha-1)))*Nstar;
    Ystar=YK*Kstar;
    Cstar=CK*Kstar;
end;

% Create instance of the DSGE class
nx=1;
nz=1;
if ds;
    ny=4;
    nu=3;
    Par=[astar;alpha;beta;delta;eta;theta];
    Names={'Capital','Growth Factor Shock','Growth Factor','Output','Consumption','Hours'};
    if log;
        v=[reallog(kstar);0;a;reallog(ystar);reallog(cstar);reallog(Nstar)];
    else
        v=[kstar;0;astar;ystar;cstar;Nstar];
    end;
else
    ny=3;
    nu=2;
    Par=[alpha;beta;delta;eta;theta];
    Names={'Capital','log of TFP','Output','Consumption','Hours'};
    if log;
        v=[reallog(Kstar);0;reallog(Ystar);reallog(Cstar);reallog(Nstar)];
    else
        v=[Kstar;0;Ystar;Cstar;Nstar];
    end;
end;

EM=DSGE(nx,ny,nz,nu,v,'Names',Names);

% Modify default settings
EM.numeric=true;
EM.reduced=true;
EM.order=2;
EM.Plot=true;
EM.Print=true;

if ds;
    EM.ds=true;    
    EM.Var(1).Xi=1;
    EM.Var(4).Xi=1;
    EM.Var(5).Xi=1;
    EM.Var(2).Plotno=1;
    EM.Var(3).Plotno=2;
    EM.Var(4).Plotno=2; EM.Var(4).Print=1; EM.Var(4).Rel=1; EM.Var(4).Corr=1;
    EM.Var(5).Plotno=2;EM.Var(5).Print=1;
    EM.Var(6).Plotno=2;EM.Var(6).Print=1;
    if log;
        EM.log=true;
        EM.outfile='EM_Numeric_ds_log';
        EM.Equations=@(x,en)EM_Sys_ds_log(x,en,Par);
    else
        EM.outfile='EM_Numeric_ds';
        EM.Equations=@(x,en)EM_Sys_ds(x,en,Par);
    end;
else
    EM.Var(2).Plotno=1; % plot shock to panel 1
    EM.Var(3).Plotno=2; EM.Var(3).Print=1; EM.Var(3).Rel=1; EM.Var(3).Corr=1;
    EM.Var(4).Plotno=2; EM.Var(4).Print=1;
    EM.Var(5).Plotno=2; EM.Var(5).Print=1;
    if log;
        EM.log=true;
        EM.outfile='EM_Numeric_log';
        EM.Equations=@(x,en)EM_Sys_log(x,en,Par);
    else
        EM.outfile='EM_Numeric';
        EM.Equations=@(x,en)EM_Sys(x,en,Par);
    end;
end;

% Transition of shocks
EM.Rho=rho;
EM.Omega=sigma;

% Solve Model
[Solution,rc]=SolveModel(EM);
    
% Simulate Model
if rc>0; 
        EM.Messages{rc}
else     
        nvar=length(EM.Var);
        for i=1:nvar
            EM.Var(i).Bound(1)=0;
            EM.Var(i).Bound(2)=5*EM.Var(i).Star;
            EM.CheckBounds=true;
        end;
        [irf,sx,rx]=SimulateModel(Solution,EM);
end;

