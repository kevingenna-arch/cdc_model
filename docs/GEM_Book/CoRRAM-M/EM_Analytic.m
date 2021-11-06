% EM_Analytic.m
%
%{
 
    Author:         Alfred Maußner
    First Version:  17 January 2017
    Revison:        02 February 2017
    Purpose:

    This Matlab script tests and solves the example model from the CoRRAM-M users guide using
    analytic differentiation.

%}

% switches
log=false; % model in logs or levels
ds=false;  % model trend stationary or difference stationary

% Calibration and stationary solution

if ~ds;
    Par=struct('Symbol','','Value',0);
    Par(1).Symbol='alpha'; Par(1).Value=0.3;
    Par(2).Symbol='beta';  Par(2).Value=0.99;
    Par(3).Symbol='delta'; Par(3).Value=0.025;
    Par(4).Symbol='eta';   Par(4).Value=2.0;
    Par(5).Symbol='theta'; Par(5).Value=2.07;

    rho=0.95;
    sigma=0.007;

    % Stationary solution
    for i=1:1:5;
        assignin('base',Par(i).Symbol,Par(i).Value);
    end;

    YK=(1-beta*(1-delta))/(alpha*beta);
    CK=YK-delta;
    temp=(theta/(1-alpha))*(CK/YK);
    Nstar=1/(1+temp);
    Kstar=(YK^(1/(alpha-1)))*Nstar;
    Ystar=YK*Kstar;
    Cstar=CK*Kstar;
else
    Par=struct('Symbol','','Value',0);
    Par(1).Symbol='alpha'; Par(1).Value=0.3;
    Par(2).Symbol='beta';  Par(2).Value=0.99;
    Par(3).Symbol='delta'; Par(3).Value=0.025;
    Par(4).Symbol='eta';   Par(4).Value=2.0;
    Par(5).Symbol='theta'; Par(5).Value=2.07;
    Par(6).Symbol='astar'; % value assigned below

    a=0.004;
    rho=0.0;
    sigma=0.018;

    % Stationary solution
    for i=1:1:5;
        assignin('base',Par(i).Symbol,Par(i).Value);
    end;
    astar=exp(a);
    Par(6).Value=astar;
    yk=(astar^eta-beta*(1-delta))/(alpha*beta);
    ck=yk-(astar-1+delta);
    temp=(theta/(1-alpha))*(ck/yk);
    Nstar=1/(1+temp);
    kstar=astar*(yk^(1/(alpha-1)))*Nstar;
    ystar=yk*kstar;
    cstar=ck*kstar;
end;

% Create instance of the DSGE class
nx=1;
nz=1;
if ~ds;    
    ny=3;
    nu=2;
    Names={'Capital','log of TFP','Output','Consumption','Hours'};
    if log;        
        Symbols={'k','z','y','c','n'};
        v=[reallog(Kstar);0;reallog(Ystar);reallog(Cstar);reallog(Nstar)];
    else        
        Symbols={'K','z','Y','C','N'};
        v=[Kstar;0;Ystar;Cstar;Nstar];
    end;
else
    ny=4;
    nu=3;
    Names={'Capital','Growth Factor Shock','Growth Factor','Output','Consumption','Hours'};
    Symbols={'k','z','a','y','c','N'};
    if log;
        v=[reallog(kstar);0;reallog(astar);reallog(ystar);reallog(cstar);reallog(Nstar)];
    else
        v=[kstar;0;astar;ystar;cstar;Nstar];
    end;
end;

EM=DSGE(nx,ny,nz,nu,v,Symbols,'Names',Names);

% Modify default settings, relevant for all cases
EM.order=3;
EM.reduced=true;
EM.Plot=true;
EM.Print=true;
EM.Rho=rho;
EM.Omega=sigma;
EM.Skew=0;

% modify default setting for the different cases
if ~ds;    
    EM.Var(2).Plotno=1; 
    EM.Var(3).Plotno=2; EM.Var(3).Print=1; EM.Var(3).Rel=1; EM.Var(3).Corr=1;
    EM.Var(4).Plotno=2; EM.Var(4).Print=1;
    EM.Var(5).Plotno=2; EM.Var(5).Print=1;
    if log;
        EM.log=true;
        EM.outfile='EM_Analytic_log';
        EM.Equations='EM_Eqs_log';
    else
        EM.outfile='EM_Analytic';
        EM.Equations='EM_Eqs';
    end;
else
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
        EM.outfile='EM_Analytic_ds_log';
        EM.Equations='EM_Eqs_ds_log';
    else
        EM.outfile='EM_Analytic_ds';
        EM.Equations='EM_Eqs_ds';
    end;
end;

% Solve Model
[Hmat,rc]=SolveModel(EM,Par);
    
% Simulate Model
if rc>0; 
        EM.Messages{rc}
else
        [irf,sx,rx]=SimulateModel(Hmat,EM);
end;

