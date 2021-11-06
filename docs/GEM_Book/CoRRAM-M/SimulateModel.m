function [irf,sx,rx] = SimulateModel(S,Model)
                                       
% Simulate DSGE model

%{
    Copyright:   Alfred Maußner

    First version:  17 March 2015
    Revisons:       29 June  2015 (employs vectorized policy functions up to third-order)
                    14 August 2015 (bug fixed, the program did not correctly stop if the size of the matrix with random numbers is incorrec).
                    26 August 2015 check boundaries added, policy functions simplified
                    24 September 2015, bug in line 336 fixed
                    25 September 2015, do not compute impulse responeses if inobs=0
                    02 October 2015, do not compute second moments if nobs=0
                    13 November 2015, ordering changed to x,z,y!
                    26 January 2016, bug in extracting results from simulated moments fixed
                    02 February 2016, bug in line 257 fixed (read bounds from Var), bug in counting invalid simulations fixed
                    12 January 2017, adapted to the DSGE class
                    18 January 2017, bug in impulse respones for ds option fixed
                    20 January 2017, HPF instead of hpfilter used to filter data
                    30 January 2017, ordering in st changed to [xt zt yt] to get correct tables
                    02 March 2018, bug in computation of relative deviations in case of zero stationary values fixed
    
    Purpose:  . Compute and display impulse responses of a DSGE model,
              . compute and print second moments of a DSGE model
           
 Input: S: structure with fields:

           . Hx_w:    nx by nw matrix, the linear part of the policy function
                      for the endogenous state variables, nw=nx+nz
           . Hy_w:    ny by nw matrix, the linear part of the policy function
                      for the control (i.e., not predetermined) variables
           . Hx_ww:   nx*nw by nw matrix the quadratic part w.r.t. states and shocks of the policy function
                      for the endogenous state variables
           . Hy_ww:   ny*nw by nw matrix, the quadratic part w.r.t. states and shocks of the policy function
                      for the control variables
           . Hw_ss:   nx by 1 matrix, the quadratic part w.r.t. the perturbation
                      parameter of the policy function for the endogenous
                      state variables.
           . Hy_ss:   ny by 1 matrix, the quadratic part w.r.t. the perturbation
                      parameter of the policy function for the control variables.
           . Hx_www:  nx*nw^2 by nw matrix, the cubic part w.r.t. states
                      and shocks of the policy function for the endogenous
                      state variables.
           . Hy_www:  ny*nw^2 by nw matrix, the cubic part w.r.t states
                      and shocks of the policy function for the control
                      variables.
           . Hw_sss:  nx by 1 matrix, the cubic part w.r.t. the perturbation parameter of the
                      policy function for the endogenous state variables
           . Hy_sss:  ny by 1 matrix, the cubic part w.r.t. the perturbation parameter of the
                      policy function for the control variables
           . Hx_ssw:  nx*nw by 1 matrix, the cubic part w.r.t. the perturbation parameter and w.r.t
                      the states and shocks of the policy function for the endogenous state variables
           . Hy_ssw:  ny*nw by 1 matrix, the cubic part w.r.t the perturbation parameter and w.r.t
                      the states and shocks of the policy function for the control variables
           . Model:   an instance of the DSGE class

%}

% open Logfile for error messages
lfh=fopen(strcat(Model.outfile,'_LogFile.txt'),'a+');

hpf=Model.hp;

% get necessary information
nx=Model.nx;
nz=Model.nz;
ny=Model.ny;
Rho=Model.Rho;
Omega=Model.Omega;
inobs=Model.inobs;  % number of periods for impulse responses
nobs=Model.nobs;    % number of periods in simulated time series
nofs=Model.nofs;    % number of simulated time series

% store stationary values 
x0=zeros(nx,1);
y0=zeros(ny,1);
w0=zeros(nx+ny,1);
m0=100*ones(inobs,nx+ny+1);

for i1=1:1:nx;
    x0(i1)=Model.Var(i1).Star;
    if x0(i1)==0; w0(i1)=1; m0(:,i1)=1; else w0(i1)=x0(i1); end
end;

for i1=1:1:ny;
    y0(i1)=Model.Var(nx+nz+i1).Star;
    if y0(i1)==0; w0(nx+i1)=1; else w0(nx+i1)=y0(i1); end
end;

% store scaling factors
if Model.ds;
    Xi=zeros(nx+ny,1);
    na=nx+ny+nz;
    i2=1;
    for i1=1:na;
        if ~strcmp(Model.Var(i1).Type,'z'); Xi(i2)=Model.Var(i1).Xi; i2=i2+1; end;
    end;
end;
% Compute and plot impulse responses
if inobs>0;    
    irf=zeros(inobs,nx+ny+1,nz);
    for i1=1:1:nz;
        irf(:,:,i1)=Impulse(i1);
        if Model.Plot; MyPlot(irf(:,:,i1),Model,i1); end;
    end;
else
    irf=0;
end;

% Compute second moments
if nobs>0;
    [sx,rx,rc,halt]=Moments(nobs,nofs);  
    if halt; fclose(lfh); return; end;

    if rc==0;
        fprintf(lfh,'No valid simulations.\n');
        fclose(lfh);
        return;
    elseif Model.Print;
        switch Model.table;
            case 'X';
                TableX(sx,rx,rc,Model);
            case 'A';
                TableA(sx,rx,rc,Model);                
            case 'B';
               TableX(sx,rx,rc,Model);
               TableA(sx,rx,rc,Model);
            otherwise;
        end; 
    end;
else
    sx=0; rx=0;
end;
fclose(lfh);
   

% PFX: the policy function for the states: x=PFX(v)
%
% Input: 
%        . v:    nx+nz+1 by 1 matrix with elements [x(t)-xstar;z(t);1]
%         
% Output:
%        . x: nx by 1 matrix

function x=PFX(v) % the policy function
  
    x=x0+S.Hx_w*v(1:nx+nz);    
    if Model.order==1; return; end;
    x=x+0.5*S.Hx_ss+0.5*kron(eye(nx),v')*S.Hx_ww*v;    
    if Model.order==2; return; end;    
    x=x+(1/6)*S.Hx_sss+0.5*kron(eye(nx),v')*S.Hx_ssw + (1/6)*kron(kron(eye(nx),v'),v')*S.Hx_www*v;
    
    return;
end

% PFY: the policy function for the costates: y=PFY(v)
%
% Input: 
%        . v:    nx+nz+1 by 1 matrix with elements [x(t)-xstar;z(t);1]
%         
% Output:
%        . y: ny by 1 matrix

function y=PFY(v) % the policy function
  
    y=y0+S.Hy_w*v(1:nx+nz);    
    if Model.order==1; return; end;
    y=y+0.5*S.Hy_ss+0.5*kron(eye(ny),v')*S.Hy_ww*v;    
    if Model.order==2; return; end;    
    y=y+(1/6)*S.Hy_sss+0.5*kron(eye(ny),v')*S.Hy_ssw + (1/6)*kron(kron(eye(ny),v'),v')*S.Hy_www*v;
    
    return;
end

% Impulse: irf=Impulse(iz)
%
% Input: . inobs: scalar, the number of periods for which the impulse responses are to be computed
%        . iz:    scalar, the index of the shock
%
% Output: . irf:  inobs by 1+nx+ny matrix, the first column stores the time profile of the shock
%                 columns 2 to nx+1 the impulse response of the state variables and columns
%                 nx+2: nx+ny+1 the impulse response of the costate variables. The ordering
%                 follows the information given in Var(i).pos.

function xt=Impulse(iz)
    
% Initialize
xt=zeros(nx,inobs+1);
yt=zeros(ny,inobs+1);
zt=zeros(nz,inobs+1);
zt(iz,2)=Omega(iz,iz);
if Model.log;
    for t=2:1:inobs;
        zt(:,t+1)=Rho*zt(:,t);
        xt(:,t+1)=S.Hx_w*[xt(:,t);zt(:,t)];
        yt(:,t)  =S.Hy_w*[xt(:,t);zt(:,t)];
    end;
    xt=[xt(:,1:inobs); yt(:,1:inobs)];
else
    xt(:,1:2)=[x0,x0];
    yt(:,1:2)=[y0,y0];
    for t=2:1:inobs;
        zt(:,t+1)=Rho*zt(:,t);
        xt(:,t+1)=x0+S.Hx_w*[(xt(:,t)-x0);zt(:,t)];
        yt(:,t)  =y0+S.Hy_w*[(xt(:,t)-x0);zt(:,t)];
    end;
    xt=[xt(:,1:inobs); yt(:,1:inobs)];
    for t=1:1:inobs;
        xt(:,t)=(xt(:,t)-[x0;y0])./w0;
    end;
end;
xt=[xt;zt(iz,1:inobs)];
xt=xt';

if Model.ds; % compute impulse responses relative to the unshocked trend path of the model   
    xtt=zeros(inobs,nx+ny);
    nxy=nx+ny;
    for j=1:1:nxy; % loop over all variables
        xi=Xi(j);
        for t=2:1:inobs;
            xtt(t,j)=xt(t,j)+xi*sum(xt(1:t-1,nx+1));
        end;
    end;
    xt=[xtt,xt(:,nx+ny+1)];
end;
xt=xt.*m0;

return;
end

% Moments [sx,rx,rc,halt]=Moments(nobs,nofs)
%
% Input:  . nobs: integer, the number of periods in each simulation
%         . nobs: integer, the number of simulations
% 
% Output: . sx nx+ny+nz by 1 matrix, standard deviations of simulated time series
%         . rx 2(nx+ny+nz) by 2(nx+ny+nz) matrix, covariances between
%           matrix current and one-period lagged variables
%         . rc, integer, the number of valid simulations
%         . halt, if true, the program returns after the user has
%           pressed the stop button. This will happen, if the program does not find
%           a file with random numbers.


% look for file and import random number if ls=1 and the file exists
function [sx,rx,rc,halt]=Moments(nobs,nofs)
        
sx=0; % vector with standard deviations
sx=zeros(2*(nx+nz+ny),2*(nx+ny+nz));
rx=0; % correlation matrix
rc=0; % number of valid simulations
halt=false; % flag to stop simulations

% see whether a file with random numbers exists and should be loaded
if Model.ls;
    ok=exist(strcat(Model.loadpath,'\eps_array.mat'),'file');
    if ok==0;
        reply=menu('Array does not exist. Create it?','Yes','Stop');
        if reply==2; halt=true; return; end;
        eps=random('Normal',0,1,nz,nobs+1,nofs);
        save(strcat(Model.loadpath,'\eps_array'),'eps');
    else
        load(strcat(Model.loadpath,'\eps_array.mat'),'eps');
        [n1,n2,n3]=size(eps);
        if n1<nz || n2<(nobs+1) || n3<nofs;
            errordlg('Wrong size of eps. The program stops.');
            halt=true;
            return;
        end;
    end;    
else
    eps=random('Normal',0,1,[nz,nobs+1,nofs]);    
end;

if Model.CheckBounds;
    % read bound for x and y
    xb=zeros(nx,2);
    for i=1:nx;
        xb(i,:)=Model.Var(i).Bound;
    end;
    yb=zeros(ny,2);
    for i=1:ny;
        yb(i,:)=Model.Var(i+nx+nz).Bound;
    end;
end;
% Simulate
ivs=0; % counts invalid simulations
for itn=1:1:nofs; % loop over the number of simulations starts here    
    % Initialize    
    err=false;
    zt=zeros(nz,nobs+1);
    xt=zeros(nx,nobs+1);
    yt=zeros(ny,nobs);
    st=zeros(nobs,nx+ny+nz);
    xt(:,1)=x0;
    zt(:,1)=Omega*eps(1:nz,1,itn);
    % simulate and filter one time series
    for t=1:1:nobs; % loop over periods starts here
        v=[xt(:,t)-x0;zt(:,t)];
        zt(:,t+1)=Rho*zt(:,t)+Omega*eps(1:nz,t,itn);
        xt(:,t+1)=PFX(v);
        if Model.CheckBounds;
            if sum(xt(:,t+1)<xb(:,1))>0; Report('x',t+1,1); ivs=ivs+1; err=true; break; end;
            if sum(xt(:,t+1)>xb(:,2))>0; Report('x',t+1,2); ivs=ivs+1; err=true; break; end;
        end;
        yt(:,t)  =PFY(v);
        if Model.CheckBounds;
            if sum(yt(:,t)<yb(:,1))>0; Report('y',t,1); ivs=ivs+1; err=true; break; end;
            if sum(yt(:,t)>yb(:,2))>0; Report('y',t,2); ivs=ivs+1; err=true; break; end;
        end;
    end;
    if ~err;
        % compute scaling factor if growth is difference stationary
        if Model.ds;
            scale=zeros(1,nobs+1);
            scale(1,1)=1;
            if ~Model.log;
                for t=1:nobs;
                    scale(1,t+1)=yt(1,t)*scale(1,t); % time path of scaling factor
                end;
            else
                for t=1:nobs;
                    scale(1,t+1)=exp(yt(1,t))*scale(1,t); % time path of scaling factor
                end;
            end;
        end;
        % prepare time series
        if Model.hp==0; % no filtering of data, compute percentage deviations from trend
            if Model.ds;
                error('You must HP filter difference stationary data.');
            else
                if ~Model.log; % data are in levels and unscaled
                    for i=1:nx;
                        st(:,i)=((xt(i,1:nobs)-x0(i))./w0(i))';
                    end;
                    for i=1:ny;
                        st(:,i+nx+nz)=((yt(i,1:nobs)-y0(i))./w0(i+nx))';
                    end
                    for i=1:nz;
                        st(:,i+nx)=zt(i,1:nobs)';
                    end
                else % data are in logs and unscaled
                    for i=1:nx;
                        st(:,i)=(xt(i,1:nobs)-x0(i))';
                    end;
                    for i=1:ny;
                        st(:,i+nx+nz)=(yt(i,1:nobs)-y0(i))';
                    end
                    for i=1:nz;
                        st(:,i+nx)=zt(i,1:nobs)';
                    end
                end
            end;
        end; % no filtering
        if Model.hp>0; % filtering demanded
            if ~Model.ds && ~Model.log; % data in levels and unscaled
                for i=1:nx;                    
                    if all(xt(i,1:nobs)>0);
                        st(:,i)=HPF(log(xt(i,1:nobs)'),hpf);
                    else
                        temp1=((xt(i,1:nobs)-x0(i))./w0(i))';
                        st(:,i)=HPF(temp1,hpf);
                    end
                end
                for i=1:ny;
                    if all(yt(i,1:nobs)>0);
                        st(:,i+nx+nz)=HPF(log(yt(i,1:nobs)'),hpf);
                    else
                        temp1=((yt(i,1:nobs)-y0(i))./w0(i+nx))';
                        st(:,i+nx+nz)=HPF(temp1,hpf);
                    end
                end
                for i=1:nz;
                    st(:,i+nx)=HPF(zt(i,1:nobs)',hpf);
                end
            end
            if ~Model.ds && Model.log; % data in logs and unscaled
                for i=1:nx;
                    st(:,i)=HPF(xt(i,1:nobs)',hpf);
                end
                for i=1:ny;
                    st(:,i+nx+nz)=HPF(yt(i,1:nobs)',hpf);
                end
                for i=1:nz;
                    st(:,i+nx)=HPF(zt(i,1:nobs)',hpf);
                end
            end
            if Model.ds && ~Model.log; % data in levels and scaled
                for i=1:nx
                    temp1=xt(i,1:nobs).*(scale(1,1:nobs).^Xi(i));
                    if sum(temp1<0);
                        i
                        error('log filtering not possible');
                    end;
                    st(:,i)=HPF(log(temp1'),hpf);
                end
                for i=1:ny
                    temp1=yt(i,1:nobs).*(scale(1,1:nobs).^Xi(i+nx));
                    if sum(temp1<0);
                        i
                        error('log filtering not possible');
                    end;
                    st(:,i+nx+nz)=HPF(log(temp1'),hpf);
                end
                for i=1:nz;
                    st(:,i+nx)=HPF(zt(i,1:nobs)',hpf);
                end
            end
            if Model.ds && Model.log; % data in logs and scaled
                for i=1:nx;
                    temp1=xt(i,1:nobs)+Xi(i)*log(scale(1,1:nobs));
                    st(:,i)=HPF(temp1',hpf);
                end
                for i=1:ny;
                    temp1=yt(i,1:nobs)+Xi(i+nx)*log(scale(1,1:nobs));
                    st(:,i+nx+nz)=HPF(temp1',hpf);
                end
                for i=1:nz;
                    st(:,i+nx)=HPF(zt(i,1:nobs)',hpf);
                end
            end
        end
        
        % Computation of second moments
        %sx=sx+sqrt(var(st));
        %rx=rx+corrcoef([st(2:nobs,:),st(1:nobs-1,:)]);
        sx1=VCX([st(2:nobs,:),st(1:nobs-1,:)],0);
        sx=sx+sx1(:,:,1);
        
    end;
end;              % loop over the number of simulations ends here

rc=nofs-ivs;
if rc>0;
    %sx=(100*sx)/rc;
    %rx=rx/rc;
    sx=sx/rc;
    [sx,rx]=GetSR(sx);
    sx=sx*100;
else
    rc=0;
    sx=0;
    rx=0;
end;
return;

% GetSR computes the correlation coefficients from the covariance matrix sx
    function [sd,rc]=GetSR(x)
        nvar=size(x,2);
        rc=zeros(nvar,nvar);
        sd=sqrt(diag(x));
        for i3=1:nvar
            for j3=1:nvar
                rc(i3,j3)=x(i3,j3)/(sd(i3)*sd(j3));
            end
        end
        return
    end
% Report(typ,lu) writes an error message to the LogFile

function Report(typ,t,lu)        
       
        if typ=='x';
            mesg='State variable i= %3i';
            for ii=1:nx;
                if lu==1;
                    if xt(ii,t)<xb(ii,1); 
                        out=[ii, itn, t, xt(ii,t), xb(ii,1)];
                        break;
                    end;
                else
                    if xt(ii,t)>xb(ii,2);
                        out=[ii, itn, t, xt(ii,t), xb(ii,2)];
                        break;
                    end;
                end;
            end; 
        else
            mesg='Costate variable i=%3i'; 
            for ii=1:ny;            
                if lu==1;
                    if yt(ii,t)<yb(ii,1);
                        out=[ii, itn, t, yt(ii,t), yb(ii,1)];
                        break;
                    end;
                else
                    if yt(ii,t)>yb(ii,2);
                        out=[ii, itn, t, yt(ii,t), yb(ii,2)];
                        break;
                    end;
                end;
            end;
        end;
        mesg=strcat(mesg,' exceeds bound at iteration %4i and time %4i ',...
            ' value= %10.6f bound= %10.6f\n');        
        fprintf(lfh,mesg,out);
        return;
    end

end

end

